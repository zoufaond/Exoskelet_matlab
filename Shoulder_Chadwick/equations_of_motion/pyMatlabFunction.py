import sympy as sp

def MatlabFunction(function,fun_name,assignto,states,constants,constants_vec):
    
    text_file = open(f"{fun_name}.m","w")
    with open(f"{fun_name}.m","w") as text_file:
        # function .. = .. ()
        header = f"function {assignto} = {fun_name}(t,y,model)" 
        print(header,file=text_file)
        
        # acces to  constants in model struct
        for var_name, var_tuple in constants.items():
            print(f'{var_name} = model.{var_name};',file=text_file)
        
        # acces to vector constants in model struct
        for var_name, var_tuple in constants_vec.items():
            for i, elem in enumerate(var_tuple):
                str_res = f'{elem} = model.{var_name}(:,{i+1});'
                print(str_res,file=text_file)
                
        for i,state_name in enumerate(states):
            print(f"{str(state_name)} = y({i+1},:);", file=text_file)
                
        sub_exprs, simplified_rhs = sp.cse(function)
        for var, expr in sub_exprs:
            print('%s = %s;' % (sp.octave_code(var),sp.octave_code(expr)),file=text_file)
        print('%s' % sp.octave_code(sp.Matrix([simplified_rhs]), assign_to = assignto),file=text_file)
