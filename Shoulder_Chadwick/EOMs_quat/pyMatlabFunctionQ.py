import sympy as sp

def MatlabFunction(function,fun_name,assignto,states,body_constants,segments,other_constants={}):
    
    text_file = open(f"{fun_name}.m","w")
    with open(f"{fun_name}.m","w") as text_file:
        # function .. = .. ()
        header = f"function {assignto} = {fun_name}(t,q,u,model)" 
        print(header,file=text_file)
        
        # acces to  constants in model struct

        for var_name, var_tuple in body_constants.items():
            n = 0
            k = 0
            if not type(var_tuple) == list:
                print(f'{var_name} = model.{var_name};',file=text_file)
            
            elif len(var_tuple) == len(segments): 
                for i, elem in enumerate(var_tuple):
                    str_res = f'{var_name}{segments[i]} = model.{var_name}{segments[i]};'
                    print(str_res,file=text_file)
            else:
                for i, elem in enumerate(var_tuple):
                    str_res = f'{elem} = model.{var_name}{segments[n]}(:,{k+1});'
                    print(str_res,file=text_file)
                    k += 1
                    if (i + 1) % 3 == 0:
                        k = 0
                        n += 1
                    else:
                        pass
        
        for var_name, var_tuple in other_constants.items():
            if not type(var_tuple) == list:
                print(f'{var_name} = model.{var_name};',file=text_file)
            else:
                for i, elem in enumerate(var_tuple):
                    str_res = f'{elem} = model.{var_name}(:,{i+1});'
                    print(str_res,file=text_file)
        
        for i,joint_states in enumerate(states):
            for j, state_name in enumerate(joint_states):
                if i == 0:
                    print(f"{str(state_name)} = q({j+1},:);", file=text_file)
                else:
                    print(f"{str(state_name)} = u({j+1},:);", file=text_file)

                
        sub_exprs, simplified_rhs = sp.cse(function)
        for var, expr in sub_exprs:
            print('%s = %s;' % (sp.octave_code(var),sp.octave_code(expr)),file=text_file)
        print('%s' % sp.octave_code(sp.Matrix([simplified_rhs]), assign_to = assignto),file=text_file)
