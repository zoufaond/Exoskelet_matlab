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
    quat_new = np.zeros([num_nodes,num_coords])
    interval_value = duration/(num_nodes - 1)
    x0 = mot_struct['coords_struct']['mot_quaternion_IC'][0,0][0]
    
    for i in range(num_coords):
        cs = sc.interpolate.CubicSpline(time,quat_coords[:,i])
        quat_new[:,i] = cs(time_new)
    trajectory = quat_new.T.flatten()
    
    return trajectory, interval_value, time_new, x0

def exp_trajectory_eul(mot_struct_name,num_nodes):
    mot_struct = sc.io.loadmat(mot_struct_name)
    time = mot_struct['coords_struct']['tout'][0,0][:,0]
    duration = time[-1]
    time_new = np.linspace(0,duration,num_nodes)
    
    eul_coords = mot_struct['coords_struct']['mot_euler_mod'][0,0]
    num_coords = np.shape(eul_coords)[1]
    eul_new = np.zeros([num_nodes,num_coords])
    interval_value = duration/(num_nodes - 1)
    x0 = mot_struct['coords_struct']['mot_euler_IC'][0,0][0]
    
    for i in range(num_coords):
        cs = sc.interpolate.CubicSpline(time,eul_coords[:,i])
        eul_new[:,i] = cs(time_new)
    trajectory = eul_new.T.flatten()
    
    return trajectory, interval_value, time_new, x0

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
    
    for i,jnt in enumerate(joints):
        if jnt != 'rev':
            traj_eul[:,i*3:(i+1)*3] = quat2eul(traj_splitted[:,i*4:(i+1)*4],jnt)
        else:
            traj_eul[:,i*3] = traj_splitted[:,i*4]
            
    traj_eul = traj_eul*180/np.pi
    
#     text_file = open(f"{file_name}","w")
#     with open(f'{file_name}',"w") as text_file:
#         print(f'Simulation',file=text_file)
#         print(f'nRows={num_nodes}',file=text_file)
#         print(f'nColumns={10+5}',file = text_file)
#         print(f'endheader',file = text_file)
#         print(f'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y',file = text_file)
        
#         for i in range(num_nodes):
#             for j in range(10):
#                 if j == 0 :
#                     print(f'{time[i]}  0.000000  0.000000  0.000000  {traj_eul[i,j]}', end = '  ',file = text_file)
#                 elif j == 9:
#                     print(f'{traj_eul[i,j]}  0.000000', file = text_file)
#                 else:
#                     print(f'{traj_eul[i,j]}', end = '  ',file = text_file)
                        
        
#     text_file.close()
    
    return traj_splitted

def quat2matfile(solution,activations,num_q,num_states,num_nodes,time,file_name):
    
    trajectories = np.zeros([num_nodes, num_q])
    for i in range(num_q):
        trajectories[:,i] = solution[(i)*num_nodes:(i+1)*num_nodes]
    input_data = input2mat(solution, num_nodes, num_states, activations, time)
    data = {
                'tout': time,
                'trajectories': trajectories,
                'inputs': input_data
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
    rot = spat.from_quat(quat)
    eul = rot.as_euler(seq)
    
    return eul

def eul2quat(eul,seq):
    rot = spat.from_euler(seq,eul)
    quat = rot.as_quat(scalar_first=True)
    
    return quat