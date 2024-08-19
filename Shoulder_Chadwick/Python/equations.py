import sympy as sp
import numpy as np
import scipy as sc
import sympy.physics.mechanics as me
from scipy.spatial.transform import Rotation as spat
import pickle


def das_trajectory(data_struct,num_nodes,duration,weight):
    time = np.linspace(0.0, duration, num=num_nodes)
    interval_value = duration/(num_nodes - 1)
    x0 = data_struct['params']['InitPosOptQuat'][0,0]['initCondQuat'].item()
    x0eul = data_struct['params']['InitPosOptQuat'][0,0]['initCondEul'].item()
    GH_motion_Eul = np.array([np.ones(num_nodes)*x0eul[6],
                          (-np.cos(time*np.pi)+1)*weight+x0eul[7],
                          np.ones(num_nodes)*x0eul[8]]).T
    GH_motion_R = spat.from_euler('YZY',GH_motion_Eul)
    GH_motion_Q = GH_motion_R.as_quat(scalar_first=True).T
    
    traj = np.zeros(13*num_nodes)
    for i in range(8):
        traj[i*num_nodes:(i+1)*num_nodes] = x0[i]
    traj[8*num_nodes:(8+4)*num_nodes] = GH_motion_Q.flatten()
    traj[(8+4)*num_nodes:(8+5)*num_nodes] = x0[12]
    
    traj_split = np.vstack(np.split(traj,13))
    d_traj = np.concatenate((np.zeros((13,1)),np.diff(traj_split)),axis=1)/interval_value
    init_guess = np.concatenate((traj,d_traj.flatten()))
    
    return traj, init_guess

def polynomials(model_struct,q,derive,model_params_struct, initCond_name):
    q_thorax0 = me.dynamicsymbols('q_thorax0')
    q_thorax1 = me.dynamicsymbols('q_thorax1')
    q_thorax2 = me.dynamicsymbols('q_thorax2')
    q_thorax3 = me.dynamicsymbols('q_thorax3')
    q_radius = (me.dynamicsymbols('q_radius'))

    q_new = q.copy()
    q_new.append(q_radius)
    q_new.insert(0,q_thorax3)
    q_new.insert(0,q_thorax2)
    q_new.insert(0,q_thorax1)
    q_new.insert(0,q_thorax0)
    
    qpol = q_new[1:4]+q_new[5:8]+q_new[9:12]+q_new[13:16]+q_new[16:18]
    nmus = len(model_struct['model_simplified_quat']['muscles'][0,0][0])
    act = me.dynamicsymbols('act_1:' + str(nmus + 1))
    
    if derive == 'symbolic':
        fmax = sp.symbols('fmax_1:' + str(nmus + 1))
        lceopt = sp.symbols('lceopt_1:' + str(nmus + 1))
        lslack = sp.symbols('lslack_1:' + str(nmus + 1))
    elif derive == 'numeric':
        fmax = []
        lceopt = []
        lslack = []
        
        for i in range(nmus):
            fmax.append(model_params_struct['params'][initCond_name][0,0]['fmax'][0,0][i,0].item())
            lceopt.append(model_params_struct['params'][initCond_name][0,0]['lceopt'][0,0][i,0].item())
            lslack.append(model_params_struct['params'][initCond_name][0,0]['lslack'][0,0][i,0].item())
            
    mus_lengths = sp.zeros(nmus,1)
    mus_forces = sp.zeros(nmus,1)
    
    for imus in range(nmus):
        muscle = model_struct['model_simplified_quat']['muscles'][0,0][0,imus]
        name = muscle['name'].item()
        npolterms = muscle['lparam_count'].item()
        polcoeff_np = muscle['lcoefs'].item()
        polcoeff = sp.Matrix(polcoeff_np)
        expon_np = muscle['lparams'].item()
        expon = sp.Matrix(expon_np)
        musdof = muscle['dof_indeces'].item()
        nmusdof = muscle['dof_count'].item()
        L = 0

        for i in range(npolterms.item()):
         # Add this term's contribution to the muscle length
            term = polcoeff[i]
            for j in range(nmusdof.item()):
                mdof = musdof[j];
                for k in range(expon_np[i,j].item()):
                    term = term * qpol[int(mdof-1)];

            L = L + term;

        mus_lengths[imus] = L
        mus_forces[imus] = muscle_force(act[imus],L,fmax[imus],lceopt[imus],lslack[imus])
        
    jacobian = -mus_lengths.jacobian(qpol).T
    FQ = jacobian * mus_forces
    TEsc = invJtrans(q_new[4:8])*sp.Matrix(FQ[3:6])
    TEac = invJtrans(q_new[8:12])*sp.Matrix(FQ[6:9])
    TEgh = invJtrans(q_new[12:16])*sp.Matrix(FQ[9:12])
    TEel = sp.Matrix(FQ[12:13])
    TE = sp.Matrix(TEsc).col_join(TEac).col_join(TEgh).col_join(TEel)
    res = TE.subs(q_new[17],0)
     
    return res, act

def create_parameters_dict(model_params_struct, initCond_name):
    data_struct = model_params_struct['params']['model'][0,0]
    data_names = data_struct.dtype.names
    # data = data_struct[data_names[1]]
    element_names = []
    element_values = []
    symbolic_list = []
    # element = data_struct['offset_scapula']
    
    for name in data_names:
        element = data_struct[name].item()
        element_length = element.shape[1]
        
        if element_length == 1:
            element_names.append(name)
            element_values.append(element.item())
            
        else:
            for k in range(element_length):
                element_names.append((name+'_'+str(k+1)))
                element_values.append(element[0,k])
                
    muscles_params_struct = model_params_struct['params'][initCond_name][0,0]
    muscles_params_names = ['fmax','lceopt','lslack']
    
    muscle_names = []
    muscle_values = []
    missing = (20,31,32,33)
    
    for name in muscles_params_names:
        muscles_param = muscles_params_struct[name].item()
        muscles_params_len = muscles_param.shape[0]
        
        for k in range(muscles_params_len):
            if k == 20:
                pass
            elif k == 31:
                pass
            elif k == 32:
                pass
            elif k == 33:
                pass
            else:
                element_names.append((name+'_'+str(k+1)))
                element_values.append(muscles_param[k,0])
                
                
#     for i in range(36):
#         act_name = 'act_'+str(i+1)
        
#         # if i == 19 or i == 20:
#         #     act_val = 0.2
#         # else:
#         act_val = 0
            
        # element_names.append(act_name)
        # element_values.append(act_val)
    
    
    symbolic_list = sp.symbols(element_names)
    dict_vals = dict(zip(symbolic_list, element_values))
        
    return dict_vals, symbolic_list, element_values
    
def save_eoms_explicit(MM,FO,q,w,name):
    # Sym2Dym,Dym2Sym = create_sym_dym_dict(q,w)
    qs = sp.symbols('q1:14')
    us = sp.symbols('u1:11')
    
    states = [qs,us]
    subs_q = {q[i]: qs[i] for i in range(len(q))}
    subs_u = {w[i]: us[i] for i in range(len(w))}
    
    MM_sym = me.msubs(MM,subs_q,subs_u)
    FO_sym = me.msubs(FO,subs_q,subs_u)
        
    with open(name,'wb') as f:
        pickle.dump([MM_sym,FO_sym,qs,us],f)
        
def save_eoms_implicit(eoms_implicit,q,w,name):
    subs_q, subs_u, subs_dq, subs_du, _, _, _, _ = create_sym_dym_dict(q,w)
    eoms_subs = me.msubs(eoms_implicit,subs_q,subs_u,subs_dq,subs_du)
    
    with open(name,'wb') as f:
        pickle.dump([eoms_subs],f)
        
def load_eoms_implicit(q,w,name):
    with open(name,'rb') as f:
        eoms_impl = pickle.load(f)
        
    eoms_impl_sym = sp.Matrix(eoms_impl)
        
    _, _, _, _, subs_qback, subs_uback, subs_dqback, subs_duback = create_sym_dym_dict(q,w)
    eoms_back = me.msubs(eoms_impl_sym,subs_qback,subs_uback,subs_dqback,subs_duback)
        
    return eoms_impl_sym
    
        
def load_eoms(name):

    with open(name,'rb') as f:
        MM,FO,q,w = pickle.load(f)
        
    return MM,FO,list(q),list(w)
    
    
def create_sym_dym_dict(q,w):
    qs = sp.symbols('q0:13')
    us = sp.symbols('u1:11')
    dqs = sp.symbols('dq0:13')
    dus = sp.symbols('du1:11')
    
    subs_q = {q[i]: qs[i] for i in range(len(q))}
    subs_u = {w[i]: us[i] for i in range(len(w))}
    subs_dq = {q[i].diff(): dqs[i] for i in range(len(q))}
    subs_du = {w[i].diff(): dus[i] for i in range(len(w))}
    
    subs_qback = swtich_dict(subs_q)
    subs_uback = swtich_dict(subs_u)
    subs_dqback = swtich_dict(subs_dq)
    subs_duback = swtich_dict(subs_du)
    
    return subs_q, subs_u, subs_dq, subs_du, subs_qback, subs_uback, subs_dqback, subs_duback
    
    for i,seg in enumerate(segment):
        if joints[i] == 'quat':
            for j in ('0','1','2','3'):
                qsym.append(sp.symbols('q' + j+ '_' + seg))
        elif joints[i] == 'rotaxis':
            qsym.append(sp.symbols('q_' + seg))
        else:
            pass
        
        if joints[i] == 'quat':
            for j in ('1','2','3'):
                wsym.append(sp.symbols('w' + j + '_' + seg))
        elif joints[i] == 'rotaxis':
            wsym.append(sp.symbols('w_'+ seg))
        else:
            pass
        
    symbolic_list = qsym+wsym
    dynamic_list = qdyn+wdyn
    sym_dym_dict = dict(zip(symbolic_list,dynamic_list))
    dym_sym_dict = {y: x for x, y in sym_dym_dict.items()}
        
        
    return sym_dym_dict, dym_sym_dict
    
    
    
def create_eoms(model_struct,model_params_struct,initCond_name, derive = 'symbolic'):
    
    # symbols
    t = sp.symbols('t')

    states = ['q','w']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    joints = ['quat','quat','quat','rotaxis','weld','weld']
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    w = []
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []
    
    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
        for i,seg in enumerate(segment):
            for j in ('1','2','3'):
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))

            if joints[i] == 'quat':
                for j in ('1','2','3'):
                    w.append(me.dynamicsymbols('w' + j + '_' + seg))
            elif joints[i] == 'rotaxis':
                w.append(me.dynamicsymbols('w_'+ seg))
            else:
                pass

            if joints[i] == 'quat':
                for j in ('0','1','2','3'):
                    q.append(me.dynamicsymbols('q' + j+ '_' + seg))
            elif joints[i] == 'rotaxis':
                q.append(me.dynamicsymbols('q_' + seg))
            else:
                pass

            mass.append(sp.symbols('mass_'+seg))
            frame.append(me.ReferenceFrame('frame_' + str(seg)))
            point_offset.append(me.Point('point_offset_' + str(seg)))
            masscenter.append(me.Point('masscenter_' + str(seg)))
            
    elif derive == 'numeric':
        data_struct = model_params_struct['params']['model'][0,0]
    
        g = data_struct['g'][0,0].item()
        c = data_struct['c'][0,0].item()
        for i,seg in enumerate(segment):
            for j in range(3):
                inertia.append(data_struct['I_' + seg][0,0][0,j])
                com.append(data_struct['com_' + seg][0,0][0,j])
                offset.append(data_struct['offset_' + seg][0,0][0,j])

            if joints[i] == 'quat':
                for j in ('1','2','3'):
                    w.append(me.dynamicsymbols('w' + j + '_' + seg))
            elif joints[i] == 'rotaxis':
                w.append(me.dynamicsymbols('w_'+ seg))
            else:
                pass

            if joints[i] == 'quat':
                for j in ('0','1','2','3'):
                    q.append(me.dynamicsymbols('q' + j+ '_' + seg))
            elif joints[i] == 'rotaxis':
                q.append(me.dynamicsymbols('q_' + seg))
            else:
                pass

            mass.append(data_struct['mass_' + seg][0,0].item())
            frame.append(me.ReferenceFrame('frame_' + str(seg)))
            point_offset.append(me.Point('point_offset_' + str(seg)))
            masscenter.append(me.Point('masscenter_' + str(seg)))

    # create eoms constraints (unit quaternions)
    constraints = sp.Matrix([[q[0]**2+q[1]**2+q[2]**2+q[3]**2-1],
                         [q[4]**2+q[5]**2+q[6]**2+q[7]**2-1],
                         [q[8]**2+q[9]**2+q[10]**2+q[11]**2-1]])
    
    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = []
        for i in range(3):
            offset_thorax.append(data_struct['offset_thorax'][0,0][0,i].item())

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient(frame_ground, 'Quaternion', q[0:4])
    frame[0].set_ang_vel(frame_ground,w[0]*frame[0].x + w[1]*frame[0].y + w[2]*frame[0].z)

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]
    kindeq = []
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # orient frame w.r.t. child's frame
        frame[i].orient(frame[i-1], 'Quaternion', q[0+i*4:4+i*4])
        frame[i].set_ang_vel(frame[i-1],w[0+i*3]*frame[i].x + w[1+i*3]*frame[i].y + w[2+i*3]*frame[i].z)

        # set masscenter points 
        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set damping in joints (c * angular_velocity)
        damping = -c*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))
        # set kinematic differential equations (q_dot = f(q,u))

    for i in range(0,3):
        kindeq.append(q[0+i*4].diff(t) - 0.5 * (-w[0+i*3]*q[1+i*4] - w[1+i*3]*q[2+i*4] - w[2+i*3]*q[3+i*4]))
        kindeq.append(q[1+i*4].diff(t) - 0.5 * (w[0+i*3]*q[0+i*4] + w[2+i*3]*q[2+i*4] - w[1+i*3]*q[3+i*4]))
        kindeq.append(q[2+i*4].diff(t) - 0.5 * (w[1+i*3]*q[0+i*4] - w[2+i*3]*q[1+i*4] + w[0+i*3]*q[3+i*4]))
        kindeq.append(q[3+i*4].diff(t) - 0.5 * (w[2+i*3]*q[0+i*4] + w[1+i*3]*q[1+i*4] - w[0+i*3]*q[2+i*4]))

    # symbols for ulna
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')
    
    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = []
        EL_rot_axis = []
        PSY_rot_axis = []
        
        for i in range(3):
            offset_humerus_rot.append(data_struct['offset_humerus_rot'][0,0][0,i].item())
            EL_rot_axis.append(data_struct['EL_rot_axis'][0,0][0,i].item())
            PSY_rot_axis.append(data_struct['PSY_rot_axis'][0,0][0,i].item())
        

    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[12])

    frame[3].set_ang_vel(ulna_rot_frame,
                          w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    kindeq.append(w[9]-q[12].diff(t))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].z,0)
    frame[4].set_ang_vel(frame[3],0)
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    frame[4].ang_vel_in(frame[3])

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[3].x + offset[1+4*3]*frame[3].y + offset[2+4*3]*frame[3].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[3])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))

    # Contact between scapula and thorax
    if derive == 'symbolic':
        contTS = sp.symbols('contTS_1:4')
        contAI = sp.symbols('contAI_1:4')
        elips_trans = sp.symbols('elips_trans_1:4')
        elips_dim = sp.symbols('elips_dim_1:4')
        k_contact_in, eps_in = sp.symbols('k_contact_in eps_in')
        k_contact_out, eps_out = sp.symbols('k_contact_out eps_out')
        second_elips_dim = sp.Symbol('second_elips_dim')
    elif derive == 'numeric':
        contTS = []
        contAI = []
        elips_trans = []
        elips_dim = []
        
        for i in range(3):
            contTS.append(data_struct['contTS'][0,0][0,i].item())
            contAI.append(data_struct['contAI'][0,0][0,i].item())
            elips_trans.append(data_struct['elips_trans'][0,0][0,i].item())
            elips_dim.append(data_struct['elips_dim'][0,0][0,i].item())
        k_contact_in = data_struct['k_contact_in'][0,0].item()
        k_contact_out = data_struct['k_contact_out'][0,0].item()
        eps_in = data_struct['eps_in'][0,0].item()
        eps_out = data_struct['eps_out'][0,0].item()
        second_elips_dim = data_struct['second_elips_dim'][0,0].item()

    # contact points 
    contact_point1 = me.Point('CP1')
    contact_point1.set_pos(point_offset[0],contTS[0]*frame[1].x+contTS[1]*frame[1].y  +contTS[2]*frame[1].z)
    # contact_point1.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
    # contact_point2.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point2.v2pt_theory(point_offset[0],frame_ground,frame[1])

    ## contact forces

    # Distances between contact points and thorax frame
    x_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.x)
    y_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.y)
    z_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.z)
    x_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.x)
    y_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.y)
    z_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.z)

    # Contact forces
    f1_in = ((x_pos1-elips_trans[0])/elips_dim[0])**2+((y_pos1-elips_trans[1])/elips_dim[1])**2+((z_pos1-elips_trans[2])/elips_dim[2])**2-1
    f1_out = ((x_pos1-elips_trans[0])/(second_elips_dim*elips_dim[0]))**2+((y_pos1-elips_trans[1])/(second_elips_dim*elips_dim[1]))**2+((z_pos1-elips_trans[2])/(second_elips_dim*elips_dim[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    F1_out = 1/2*(f1_out+sp.sqrt(f1_out**2+eps_out**2))
    Fx1 = -(k_contact_in*F1_in+k_contact_out*F1_out)*(x_pos1-elips_trans[0])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[0]**2)
    Fy1 = -(k_contact_in*F1_in+k_contact_out*F1_out)*(y_pos1-elips_trans[1])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[1]**2)
    Fz1 = -(k_contact_in*F1_in+k_contact_out*F1_out)*(z_pos1-elips_trans[2])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/elips_dim[0])**2+((y_pos2-elips_trans[1])/elips_dim[1])**2+((z_pos2-elips_trans[2])/elips_dim[2])**2-1
    f2_out = ((x_pos2-elips_trans[0])/(second_elips_dim*elips_dim[0]))**2+((y_pos2-elips_trans[1])/(second_elips_dim*elips_dim[1]))**2+((z_pos2-elips_trans[2])/(second_elips_dim*elips_dim[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
    F2_out = 1/2*(f2_out+sp.sqrt(f2_out**2+eps_out**2))
    Fx2 = -(k_contact_in*F2_in+k_contact_out*F2_out)*(x_pos2-elips_trans[0])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[0]**2)
    Fy2 = -(k_contact_in*F2_in+k_contact_out*F2_out)*(y_pos2-elips_trans[1])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[1]**2)
    Fz2 = -(k_contact_in*F2_in+k_contact_out*F2_out)*(z_pos2-elips_trans[2])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[2]**2)

    # applying contact forces to contact points in thorax frame
    cont_force1 = [(contact_point1,frame_ground.x*Fx1+frame_ground.y*Fy1+frame_ground.z*Fz1)]
    cont_force2 = [(contact_point2,frame_ground.x*Fx2+frame_ground.y*Fy2+frame_ground.z*Fz2)]
    CONT = cont_force1+cont_force2

    TE,activations = polynomials(model_struct,q,derive,model_params_struct,initCond_name)
    Torques = me.dynamicsymbols('Torques1:13')
    torque = []
    for i in range(4):
        torque.append((frame[i],(Torques[i*3]*frame[i].x+Torques[i*3+1]*frame[i].y+Torques[i*3+2]*frame[i].z)*2))
    

    KM = me.KanesMethod(frame_ground, q_ind=q, u_ind=w, kd_eqs=kindeq)
    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()

    # return q,w,fr,frstar,kindeq, Torques
    return MM,FO,TE,q,w,fr,frstar,kindeq,xdot,constraints,activations



def create_eoms_u0state(model_struct,model_params_struct,initCond_name, derive = 'symbolic'):
    # symbols
    t = sp.symbols('t')

    states = ['q','w','u0']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    joints = ['quat','quat','quat','rotaxis','weld','weld']
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    w = []
    u0 = []
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []

    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
    elif derive == 'numeric':
        data_struct = model_params_struct['params']['model'][0,0]
        g = data_struct['g'][0,0].item()
        c = data_struct['c'][0,0].item()


    for i,seg in enumerate(segment):
        for idat,j in enumerate(('1','2','3')):
            if derive == 'symbolic':
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))
            elif derive == 'numeric':
                inertia.append(data_struct['I_' + seg][0,0][0,idat])
                com.append(data_struct['com_' + seg][0,0][0,idat])
                offset.append(data_struct['offset_' + seg][0,0][0,idat])

        if joints[i] == 'quat':
            for j in ('1','2','3'):
                w.append(me.dynamicsymbols('w' + j + '_' + seg))
            for j in ('0','1','2','3'):
                q.append(me.dynamicsymbols('q' + j+ '_' + seg))
            u0.append(me.dynamicsymbols('u0'+ '_' + seg))
        elif joints[i] == 'rotaxis':
            w.append(me.dynamicsymbols('w_'+ seg))
            q.append(me.dynamicsymbols('q_' + seg))
        else:
            pass

        if derive == 'symbolic':
            mass.append(sp.symbols('mass_'+seg))
        elif derive == 'numeric':
            mass.append(data_struct['mass_' + seg][0,0].item())

        frame.append(me.ReferenceFrame('frame_' + str(seg)))
        point_offset.append(me.Point('point_offset_' + str(seg)))
        masscenter.append(me.Point('masscenter_' + str(seg)))

    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = []
        for i in range(3):
            offset_thorax.append(data_struct['offset_thorax'][0,0][0,i].item())

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient(frame_ground, 'Quaternion', q[0:4])
    N_w_A = (frame[0].ang_vel_in(frame_ground))
    kinematical = sp.Matrix([
        u0[0] - q[0].diff(t),
        w[0] - N_w_A.dot(frame[0].x),
        w[1] - N_w_A.dot(frame[0].y),
        w[2] - N_w_A.dot(frame[0].z),
    ])

    for i in range(1,3):
        frame[i].orient(frame[i-1], 'Quaternion', q[0+i*4:4+i*4])
        N_w_A = (frame[i].ang_vel_in(frame[i-1]))
        kinematical = kinematical.col_join(sp.Matrix([
                            u0[i] - q[i*4].diff(t),
                             w[i*3] - N_w_A.dot(frame[i].x),
                             w[i*3+1] - N_w_A.dot(frame[i].y),
                             w[i*3+2] - N_w_A.dot(frame[i].z),
                                ]))

    frame[0].set_ang_vel(frame_ground,w[0]*frame[0].x + w[1]*frame[0].y + w[2]*frame[0].z)

    for i in range(1,3):
        frame[i].set_ang_vel(frame[i-1],w[0+i*3]*frame[i].x + w[1+i*3]*frame[i].y + w[2+i*3]*frame[i].z)

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # set masscenter points 
        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set damping in joints (c * angular_velocity)
        damping = -c*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))
        # set kinematic differential equations (q_dot = f(q,u))

    # symbols for ulna
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')

    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = []
        EL_rot_axis = []
        PSY_rot_axis = []

        for i in range(3):
            offset_humerus_rot.append(data_struct['offset_humerus_rot'][0,0][0,i].item())
            EL_rot_axis.append(data_struct['EL_rot_axis'][0,0][0,i].item())
            PSY_rot_axis.append(data_struct['PSY_rot_axis'][0,0][0,i].item())


    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[12])

    frame[3].set_ang_vel(ulna_rot_frame,
                          w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    kinematical = kinematical.col_join(sp.Matrix([w[9]-q[12].diff(t)]))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].z,0)
    frame[4].set_ang_vel(frame[3],0)
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    frame[4].ang_vel_in(frame[3])

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[3].x + offset[1+4*3]*frame[3].y + offset[2+4*3]*frame[3].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[3])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))

    # Contact between scapula and thorax
    if derive == 'symbolic':
        contTS = sp.symbols('contTS_1:4')
        contAI = sp.symbols('contAI_1:4')
        elips_trans = sp.symbols('elips_trans_1:4')
        elips_dim = sp.symbols('elips_dim_1:4')
        k_contact_in, eps_in = sp.symbols('k_contact_in eps_in')
        k_contact_out, eps_out = sp.symbols('k_contact_out eps_out')
        second_elips_dim = sp.Symbol('second_elips_dim')
    elif derive == 'numeric':
        contTS = []
        contAI = []
        elips_trans = []
        elips_dim = []

        for i in range(3):
            contTS.append(data_struct['contTS'][0,0][0,i].item())
            contAI.append(data_struct['contAI'][0,0][0,i].item())
            elips_trans.append(data_struct['elips_trans'][0,0][0,i].item())
            elips_dim.append(data_struct['elips_dim'][0,0][0,i].item())
        k_contact_in = data_struct['k_contact_in'][0,0].item()
        k_contact_out = data_struct['k_contact_out'][0,0].item()
        eps_in = data_struct['eps_in'][0,0].item()
        eps_out = data_struct['eps_out'][0,0].item()
        second_elips_dim = data_struct['second_elips_dim'][0,0].item()

    # contact points 
    contact_point1 = me.Point('CP1')
    contact_point1.set_pos(point_offset[0],contTS[0]*frame[1].x+contTS[1]*frame[1].y  +contTS[2]*frame[1].z)
    # contact_point1.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
    # contact_point2.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point2.v2pt_theory(point_offset[0],frame_ground,frame[1])

    ## contact forces

    # Distances between contact points and thorax frame
    x_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.x)
    y_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.y)
    z_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.z)
    x_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.x)
    y_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.y)
    z_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.z)

    # Contact forces
    f1_in = ((x_pos1-elips_trans[0])/elips_dim[0])**2+((y_pos1-elips_trans[1])/elips_dim[1])**2+((z_pos1-elips_trans[2])/elips_dim[2])**2-1
    f1_out = ((x_pos1-elips_trans[0])/(second_elips_dim*elips_dim[0]))**2+((y_pos1-elips_trans[1])/(second_elips_dim*elips_dim[1]))**2+((z_pos1-elips_trans[2])/(second_elips_dim*elips_dim[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    F1_out = 1/2*(f1_out+sp.sqrt(f1_out**2+eps_out**2))
    Fx1 = -(k_contact_in*F1_in+k_contact_out*F1_out)*(x_pos1-elips_trans[0])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[0]**2)
    Fy1 = -(k_contact_in*F1_in+k_contact_out*F1_out)*(y_pos1-elips_trans[1])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[1]**2)
    Fz1 = -(k_contact_in*F1_in+k_contact_out*F1_out)*(z_pos1-elips_trans[2])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/elips_dim[0])**2+((y_pos2-elips_trans[1])/elips_dim[1])**2+((z_pos2-elips_trans[2])/elips_dim[2])**2-1
    f2_out = ((x_pos2-elips_trans[0])/(second_elips_dim*elips_dim[0]))**2+((y_pos2-elips_trans[1])/(second_elips_dim*elips_dim[1]))**2+((z_pos2-elips_trans[2])/(second_elips_dim*elips_dim[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
    F2_out = 1/2*(f2_out+sp.sqrt(f2_out**2+eps_out**2))
    Fx2 = -(k_contact_in*F2_in+k_contact_out*F2_out)*(x_pos2-elips_trans[0])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[0]**2)
    Fy2 = -(k_contact_in*F2_in+k_contact_out*F2_out)*(y_pos2-elips_trans[1])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[1]**2)
    Fz2 = -(k_contact_in*F2_in+k_contact_out*F2_out)*(z_pos2-elips_trans[2])*(elips_dim[0]**2+elips_dim[1]**2+elips_dim[2]**2)/(elips_dim[2]**2)

    # applying contact forces to contact points in thorax frame
    cont_force1 = [(contact_point1,frame_ground.x*Fx1+frame_ground.y*Fy1+frame_ground.z*Fz1)]
    cont_force2 = [(contact_point2,frame_ground.x*Fx2+frame_ground.y*Fy2+frame_ground.z*Fz2)]
    CONT = cont_force1+cont_force2

    TE,activations = polynomials(model_struct,q,derive,model_params_struct,initCond_name)

    q_dep = [q[0],q[4],q[8]]
    q_ind = q[1:4]+q[5:8]+q[9:]
    holonomic = sp.Matrix([[q[0]**2+q[1]**2+q[2]**2+q[3]**2-1],
                             [q[4]**2+q[5]**2+q[6]**2+q[7]**2-1],
                             [q[8]**2+q[9]**2+q[10]**2+q[11]**2-1]])
    KM = me.KanesMethod(frame_ground, q_ind=q_ind, u_ind=w, kd_eqs=kinematical,
                        q_dependent = q_dep, u_dependent = u0, 
                        configuration_constraints=holonomic,
                        velocity_constraints=holonomic.diff(t))

    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()
    
    return MM,FO,TE,q,w,u0,fr,frstar,kinematical,xdot,holonomic,activations



def muscle_force(act, lmt, fmax, lceopt, lslack):
    lm = lmt - lslack
    f_gauss = 0.25
    kpe = 5
    epsm0 = 0.6
    fpe = (sp.exp(kpe*(lm / lceopt - 1)/epsm0)-1)/(sp.exp(kpe)-1)
    flce = (sp.exp(-(lm / lceopt - 1)**2 / f_gauss))
    force = (flce * (act*0.5) +  fpe) * fmax
    
    return force

def invJtrans(quat):
    q1 = quat[0];
    q2 = quat[1];
    q3 = quat[2];
    q4 = quat[3];
    res = sp.Matrix([[ q1/2,  q4/2, -q3/2],
                     [-q4/2,  q1/2,  q2/2],
                     [ q3/2, -q2/2,  q1/2]])

    return res

def G(quat):
    Q0 = quat[0];
    Q1 = quat[1];
    Q2 = quat[2];
    Q3 = quat[3];
    res = np.array([[-Q1, Q0, Q3, -Q2],
                    [-Q2,-Q3, Q0, Q1],
                    [-Q3, Q2, -Q1, Q0]])
    return res

def switch_dict(my_dict):
    my_new_dict = dict(zip(my_dict.values(), my_dict.keys()))
    return my_new_dict