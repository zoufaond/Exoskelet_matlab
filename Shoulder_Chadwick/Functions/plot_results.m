function res = plot_results(euler_solution_name,quaternion_solution_name,inverse_kinematics_name)

    quat_sol = load(quaternion_solution_name);
    eul_sol = load(euler_solution_name);
    ik_sol = load(inverse_kinematics_name);

end