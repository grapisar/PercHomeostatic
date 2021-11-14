from sys import argv

def input_params():
    check = True
    try:
        ins = argv[1]
        dist = argv[2]
        SAMPLES = argv[3]
        outstring = argv[4]
        plot = argv[5]
    except IndexError:
        print("Missing inputs. Please insert all required inputs. See README.md for details.")
        check = False
    else:
        if (dist != 'gauss' and dist != 'uniform'):
            print("Invalid weight distribution. Please insert 'uniform' or 'gauss'.")
            check = False
        if (plot != 'plot' and plot != 'noplot'):
            print("Invalid plot argument. Insert 'plot' to visualize results or 'noplot' to avoid visualization.")
            check = False
        try:
            SAMPLES = int(SAMPLES)
            open(ins)
        except ValueError as e1:
            print(e1)
            check = False
        except FileNotFoundError as e2:
            print(e2)
            check = False

    return [ins,dist,SAMPLES,outstring,plot,check]