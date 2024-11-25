import sympy as sp
import numpy as np
import scipy as sc
import sympy.physics.mechanics as me
from scipy.spatial.transform import Rotation as spat
import pickle

def exp_trajectory_quat(mot_struct_name,num_nodes):
    mot_struct = sc.io.loadmat(mot_struct_name)
    time = mot_struct['coords_struct']['tout'][0,0][:,0]
    duration = time[-1]
    time_new = np.linspace(0,duration,num_nodes)
    
    quat_coords = mot_struct['coords_struct']['mot_quaternion_mod'][0,0]
    num_coords = np.shape(quat_coords)[1]
    trajectory = np.zeros([num_coords,num_nodes])
    interval_value = duration/(num_nodes - 1)
    x0 = mot_struct['coords_struct']['mot_quaternion_IC'][0,0][0]
    
    for i in range(num_coords):
        cs = sc.interpolate.CubicSpline(time,quat_coords[:,i])
        trajectory[i,:] = cs(time_new)
    
    return trajectory, interval_value, time_new, x0

def exp_trajectory_quat_2_myobj(trajectory):
    new_traj = trajectory.copy()
    num_nodes = np.shape(trajectory)[1]

    for i in range(num_nodes):
        new_traj[4:8,i] = mulQuat_np(trajectory[0:4,i],trajectory[4:8,i]).reshape(4,)

    return new_traj

def exp_trajectory_eul(mot_struct_name,num_nodes):
    mot_struct = sc.io.loadmat(mot_struct_name)
    time = mot_struct['mot_struct']['time'][0,0][:,0]
    duration = time[-1]
    time_new = np.linspace(0,duration,num_nodes)
    
    eul_coords = mot_struct['mot_struct']['euler'][0,0] #mot_euler_mod
    num_coords = np.shape(eul_coords)[1]
    eul_new = np.zeros([num_nodes,num_coords])
    interval_value = duration/(num_nodes - 1)
    # x0 = mot_struct['mot_struct']['mot_euler_IC'][0,0][0]
    
    for i in range(num_coords):
        cs = sc.interpolate.CubicSpline(time,eul_coords[:,i])
        eul_new[:,i] = cs(time_new)
    trajectory = eul_new.T.flatten()
    
    return trajectory, interval_value, time_new #, x0

def das_trajectory(data_struct,num_nodes,duration,weight, coords):
    time = np.linspace(0.0, duration, num=num_nodes)
    interval_value = duration/(num_nodes - 1)
    x0 = data_struct['params']['InitPosOptQuat'][0,0]['initCondQuat'].item()
    x0eul = data_struct['params']['InitPosOptQuat'][0,0]['initCondEul'].item()
    GH_motion_Eul = np.array([np.ones(num_nodes)*x0eul[6],
                          (-np.cos(time*np.pi)+1)*weight+x0eul[7],
                          np.ones(num_nodes)*x0eul[8]]).T
    GH_motion_R = spat.from_euler('YZY',GH_motion_Eul)
    GH_motion_Q = GH_motion_R.as_quat(scalar_first=True).T
    
    
    
    if coords == 'quaternion':
        traj = np.zeros(13*num_nodes)
        for i in range(8):
            traj[i*num_nodes:(i+1)*num_nodes] = x0[i]

            traj[8*num_nodes:(8+4)*num_nodes] = GH_motion_Q.flatten()
            traj[(8+4)*num_nodes:(8+5)*num_nodes] = x0[12]

            traj_split = np.vstack(np.split(traj,13))
            d_traj = np.concatenate((np.zeros((13,1)),np.diff(traj_split)),axis=1)/interval_value
            init_guess = np.concatenate((traj,d_traj.flatten()))
        
    elif coords == 'euler':
        traj = np.zeros(10*num_nodes)
        for i in range(6):
            traj[i*num_nodes:(i+1)*num_nodes] = x0eul[i]

            traj[6*num_nodes:(6+3)*num_nodes] = GH_motion_Eul.T.flatten()
            traj[(6+3)*num_nodes:(6+4)*num_nodes] = x0eul[9]

            traj_split = np.vstack(np.split(traj,10))
            d_traj = np.concatenate((np.zeros((10,1)),np.diff(traj_split)),axis=1)/interval_value
            init_guess = np.concatenate((traj,d_traj.flatten()))
    
    return traj, init_guess

def sol2mot_quat(solution, num_nodes, num_q, time, file_name = 'traj_opt.mot'):
    traj_quat = solution[:(num_nodes*num_q)]
    traj_splitted = np.vstack(np.split(traj_quat,num_q)).T
    joints = ('YZX','YZX','YZY','rev')
    traj_eul = np.zeros([num_nodes,10])

    # print(np.shape(traj_splitted))
    
    for i,jnt in enumerate(joints):
        if jnt != 'rev':
            for j in range(num_nodes):
                traj_eul[j,i*3:(i+1)*3] = quat2eul(traj_splitted[j,i*4:(i+1)*4],jnt)
        else:
            traj_eul[:,i*3] = traj_splitted[:,i*4]
            
    traj_eul = traj_eul*180/np.pi
    
    text_file = open(f"{file_name}","w")
    with open(f'{file_name}',"w") as text_file:
        print(f'Simulation',file=text_file)
        print(f'nRows={num_nodes}',file=text_file)
        print(f'nColumns={10+5}',file = text_file)
        print(f'endheader',file = text_file)
        print(f'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y',file = text_file)
        
        for i in range(num_nodes):
            for j in range(10):
                if j == 0 :
                    print(f'{time[i]}  0.000000  0.000000  0.000000  {traj_eul[i,j]}', end = '  ',file = text_file)
                elif j == 9:
                    print(f'{traj_eul[i,j]}  0.000000', file = text_file)
                else:
                    print(f'{traj_eul[i,j]}', end = '  ',file = text_file)
                        
        
    text_file.close()

    print('Saved to .mot file')

def sol2struct(solution,activations,num_q,num_states,num_nodes,time,num_iter_sol,time2sol,file_name):
    
    trajectories = np.zeros([num_nodes, num_q])
    for i in range(num_q):
        trajectories[:,i] = solution[(i)*num_nodes:(i+1)*num_nodes]
    inputs = input2mat(solution, num_nodes, num_states, activations, time)
    data = {
                'tout': time,
                'trajectories': trajectories,
                'inputs': inputs,
                'time2sol': time2sol,
                'num_iter_sol': num_iter_sol
                    }
    sc.io.savemat(f'{file_name}', {'data': data})
    print('Saved to .mat file')

def input2mat(solution, num_nodes, num_states, activations, time):
    activations_str = [str(activations[x]).replace('(t)','') for x in range(len(activations))]
    act_index = np.linspace(0,137,138,dtype=int)
    dict_act_index = dict(zip(sorted(activations_str), act_index))
    
    data = np.zeros([num_nodes, 137])
    
    for i in range(137):
        try:
            ind = dict_act_index.get(f'act_{i+1}')
            data[:,i] = solution[(num_states+ind)*num_nodes:(num_states+ind+1)*num_nodes]
        except:
            continue
    # timeseries_data = {
    #             'time': time,
    #             'data': data
    #                 }
    # sc.io.savemat('../timeseries_data.mat', {'timeseries_data': timeseries_data})
    # print('Saved to .mat file')
    return data

def quat2eul(quat,seq):
    rotm = Qrm(quat)

    if seq == 'YZY':
        z = np.arccos(rotm[1,1])
        yy = np.arctan2(rotm[1,2],rotm[1,0])
        y = np.arctan2(rotm[2,1],-rotm[0,1])
        eul = np.array([y,z,yy])
    elif seq == 'YZX':
        z = np.arcsin(rotm[1,0])
        x = np.arctan2(-rotm[1,2],rotm[1,1])
        y = np.arctan2(-rotm[2,0],rotm[0,0])
        eul = np.array([y,z,x])

    return eul

def eul2quat(eul,seq):
    rot = spat.from_euler(seq,eul)
    quat = rot.as_quat(scalar_first=True)
    
    return quat

def Qrm(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  np.array([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0             ,1]])

    return res

def mulQuat_sp(qa,qb):
    res = sp.Matrix([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res

def mulQuat_np(qa,qb):
    res = np.array([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res


def Qrm(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  sp.Matrix([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0             ,1]])

    return res

def YZY_seq(phi1,phi2,phi3):
    res = R_y(phi1) * R_z(phi2) * R_y(phi3)

    return res

def T_x(x):
    trans_x = sp.Matrix([[1,0,0,x],
                         [0,1,0,0],
                         [0,0,1,0],
                         [0,0,0,1]])
    return trans_x

def T_y(y):
    trans_y = sp.Matrix([[1,0,0,0],
                         [0,1,0,y],
                         [0,0,1,0],
                         [0,0,0,1]])
    return trans_y

def T_y(z):
    trans_y = sp.Matrix([[1,0,0,0],
                         [0,1,0,0],
                         [0,0,1,z],
                         [0,0,0,1]])
    return trans_y

def R_x(phix):
    rot_phix = sp.Matrix([[1,0          ,0           ,0],
                         [0,sp.cos(phix),-sp.sin(phix),0],
                         [0,sp.sin(phix), sp.cos(phix),0],
                         [0,0          ,0           ,1]])
    return rot_phix

def R_y(phiy):
    rot_phiy = sp.Matrix([[sp.cos(phiy) ,0,sp.sin(phiy),0],
                         [0           ,1,0          ,0],
                         [-sp.sin(phiy),0,sp.cos(phiy),0],
                         [0           ,0,0           ,1]])
    return rot_phiy

def R_z(phiz):
    rot_phiz = sp.Matrix([[sp.cos(phiz),-sp.sin(phiz),0,0],
                        [sp.sin(phiz), sp.cos(phiz),0,0],
                        [0           ,0            ,1,0],
                        [0           ,0            ,0,1]])
    return rot_phiz

def position(x,y,z):
    r = sp.Matrix([x,y,z,1])
    return r

def G(quat):
    Q0 = quat[0];
    Q1 = quat[1];
    Q2 = quat[2];
    Q3 = quat[3];
    res = np.array([[-Q1, Q0, Q3, -Q2],
                    [-Q2,-Q3, Q0, Q1],
                    [-Q3, Q2, -Q1, Q0]])
    return res