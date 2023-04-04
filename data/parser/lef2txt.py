# Program to translate lef files into txt with ICCAD2022 contest format
import pdb
import sys, getopt
# site
# Nangate15: 128 x 1536
# Nangate45: 380 x 2800

# TODO: for PI PO MC index

ENABLE_IO = False
ENABLE_REPLACE = True
FORCE_NAME = True

cell_infos = {} # pin_infos, pin_num
inst_infos = {}
net_infos = {}


def parseLef(file_name:str, tech_name:str): # tech_name is TA TB
    cell_flag = False
    pin_flag = False
    MC_cnt = 1
    P_cnt = 1

    cell_infos[tech_name] = {}
    with open(file_name, 'r') as f_in:
        line = f_in.readline()
        cnt = 1
        while line:
            print("Line {}: {}".format(cnt, line.strip()))
            lines = line.split()
            if(cell_flag):
                
                if ("SIZE" in line):
                    x = lines[1]
                    y = lines[3]
                    if(tech_name == "TA"):
                        cell_infos[tech_name][cell_name] = {"size_x": int(float(x) * 2000), 
                                                            "size_y": int(float(y) * 2000),
                                                            "pins":{},
                                                            "alias": "MC"+str(MC_cnt)}
                    elif(tech_name == "TB"):
                        if((cell_name in cell_infos["TA"].keys())):
                            alias = cell_infos["TA"][cell_name]["alias"]
                        else:
                            alias = ""
                        cell_infos[tech_name][cell_name] = {"size_x": int(float(x) * 2000), 
                                                            "size_y": int(float(y) * 2000),
                                                            "pins":{},
                                                            "alias": alias}
                if ("PIN" in line):
                    pin_flag = True
                    pin_name = lines[1]
                
                if ("DIRECTION" in line and pin_flag):
                    # donot consider INOUT pin
                    if ("INPUT" not in line and "OUTPUT" not in line):
                        pin_flag = False
                        pin_name = ""
                
                if(pin_flag and " RECT " in line):
                    # TODO: take multiple rect into consideration
                    x_l = float(lines[1])
                    y_l = float(lines[2])
                    x_r = float(lines[3])
                    y_r = float(lines[4])
                    pin_x = (x_l + x_r) / 2.0
                    pin_y = (y_l + y_r) / 2.0
                    if(tech_name == "TA"):
                        cell_infos[tech_name][cell_name]["pins"][pin_name] = {"loc_x": int(pin_x * 2000), "loc_y": int(pin_y * 2000), "alias": "P"+str(P_cnt)}
                    elif(tech_name == "TB"):
                        if((cell_name in cell_infos["TA"].keys())):
                            alias = cell_infos["TA"][cell_name]["pins"][pin_name]["alias"]
                        else:
                            alias = ""
                        cell_infos[tech_name][cell_name]["pins"][pin_name] = {"loc_x": int(pin_x * 2000), "loc_y": int(pin_y * 2000), "alias": alias}
                    pin_flag = False
                    P_cnt += 1

                if(cell_name in line and "END" in line):
                    if(len(cell_infos[tech_name][cell_name]["pins"]) < 2):
                        del cell_infos[tech_name][cell_name]
                        MC_cnt -= 1
                    cell_flag = False
                    MC_cnt += 1

            else:
                if ("MACRO" in line):
                    cell_flag = True
                    cell_name = line.split()[1]
                    P_cnt = 1
            
            line = f_in.readline()
            cnt += 1

def parseDef(file_name:str):
    comp_flag = False
    pins_flag = False
    nets_flag = False
    net_flag = False

    C_cnt = 1

    with open(file_name, 'r') as f_in:
        line = f_in.readline()
        cnt = 1
        while line:
            print("Line {}: {}".format(cnt, line.strip()))

            if("DIEAREA" in line):
                lines = line.split()
                d_x = float(lines[6])
                d_y = float(lines[7])
            elif("COMPONENTS" in line):
                comp_flag = "END" not in line
            elif("NETS" in line):
                nets_flag = "END" not in line
            elif("PINS" in line):
                pins_flag = "END" not in line
            elif(comp_flag):
                if("-" in line):
                    lines = line.split()
                    if(ENABLE_REPLACE):
                        lines[1] = lines[1].replace('/','_')
                        lines[1] = lines[1].replace('\[', '_')
                        lines[1] = lines[1].replace('\]', '_')
                    inst_infos[lines[1]] = {"type":lines[2],
                                            "alias":"C"+str(C_cnt)}
                    C_cnt += 1

            elif(ENABLE_IO):
                if(pins_flag):
                    if("-" in line):
                        lines = line.split()
                        if("INPUT" in line):
                            inst_infos[lines[1]] = {"type":"PI", 
                                                    "alias":"C"+str(C_cnt)}
                            C_cnt += 1

                        elif("OUTPUT" in line):
                            inst_infos[lines[1]] = {"type":"PO",
                                                    "alias":"C"+str(C_cnt)}
                            C_cnt += 1
    
            elif(nets_flag):
                if("-" in line):
                    lines = line.split()
                    net_flag = True
                    if(ENABLE_REPLACE):
                        lines[1] = lines[1].replace('/', '_')
                        lines[1] = lines[1].replace('\[', '_')
                        lines[1] = lines[1].replace('\]', '_')
                    net_name = lines[1]
                    pin_vec = []
                elif(";" in line):
                    net_flag = False
                    if(len(pin_vec) != 0):
                        if(ENABLE_REPLACE):
                            net_name = net_name.replace('/', '_')
                            net_name = net_name.replace('[', '_')
                            net_name = net_name.replace(']', '_')
                        net_infos[net_name] = pin_vec
                    net_name = ""
                    pin_vec = []
                elif(net_flag):
                    lines = line.split(") (")
                    lines[0] = lines[0].split("(")[1]
                    lines[-1] = lines[-1].split(")")[0]
                    for p_info in lines:
                        p_infos = p_info.split()
                        assert(len(p_infos) == 2)
                        if(p_infos[0] == "PIN"):
                            if(ENABLE_IO):
                                pin_vec.append(p_infos[1])
                            else:
                                continue
                        else:
                            if(ENABLE_REPLACE):
                                p_infos[0] = p_infos[0].replace('/', '_')
                                p_infos[0] = p_infos[0].replace('\[', '_')
                                p_infos[0] = p_infos[0].replace('\]', '_')
                            pin_vec.append(p_infos[0]+'/'+p_infos[1])

            line = f_in.readline()
            cnt += 1
        return d_x, d_y

def diesize(x, y):
    lcm_y = 268800
    target_area = 0.7 * x * y # 2 / sqrt(5)
    new_x = 0
    new_y = 0

    for y_mul in range(1,100):
        new_y = y_mul * lcm_y
        new_x = target_area / new_y
        if(abs(new_x - new_y)/new_y < 0.8 and abs(new_x - new_y)/new_x < 0.8):
            break
    
    d_x = int(new_x)
    d_y = int(new_y)

    return d_x, d_y

def output(output_file:str, d_x, d_y):
    f = open(output_file, "w")

    tech_num = len(cell_infos)
    f.write("NumTechnologies {}\n".format(tech_num))
    for tech, infos in sorted(cell_infos.items()):
        if(ENABLE_IO):
            f.write("Tech {} {}\n".format(tech, len(infos)+2))
        else:
            f.write("Tech {} {}\n".format(tech, len(infos)))
        for cell, cell_info in sorted(infos.items()):
            f.write("LibCell {} {} {} {}\n".format(cell_info["alias"], cell_info["size_x"], cell_info["size_y"], len(cell_info["pins"])))
            for pin, pin_info in sorted(cell_info["pins"].items()):
                f.write("Pin {} {} {}\n".format(pin_info["alias"], pin_info["loc_x"], pin_info["loc_y"]))
        if(ENABLE_IO):
            f.write("LibCell MC{} 1 1 1\n".format(len(infos)+1)) # PI
            f.write("Pin P1 0 0\n")
            f.write("LibCell MC{} 1 1 1\n".format(len(infos)+2)) # PO
            f.write("Pin P1 0 0\n")

    f.write("\n")

    f.write("DieSize 0 0 {} {}\n".format(d_x, d_y))

    f.write("\n")

    f.write("TopDieMaxUtil 80\n")
    f.write("BottomDieMaxUtil 85\n")

    f.write("\n")

    cnt_top = int(d_y / 1536)
    cnt_bottom = int(d_y / 2800)
    
    f.write("TopDieRows 0 0 {} 1536 {}\n".format(d_x, cnt_top))
    f.write("BottomDieRows 0 0 {} 2800 {}\n".format(d_x, cnt_bottom))

    f.write("\n")

    f.write("TopDieTech TA\n")
    f.write("BottomDieTech TB\n")

    f.write("\n")

    f.write("TerminalSize 1000 1000\n")
    f.write("TerminalSpacing 800\n")

    f.write("\n")

    f.write("NumInstances {}\n".format(len(inst_infos)))
    for inst_name, inst_info in inst_infos.items():
        f.write("Inst {} {}\n".format(inst_info["alias"], cell_infos["TA"][inst_info["type"]]["alias"]))

    f.write("\n")

    f.write("NumNets {}\n".format(len(net_infos) - 1))
    net_cnt = 1
    clk_flag = True
    for net_name, pin_infos in net_infos.items():
        if(clk_flag):
            clk_flag = False
            continue
        f.write("Net {} {}\n".format("N"+str(net_cnt), len(pin_infos)))
        net_cnt += 1
        for p_name in pin_infos:
            inst_name, pin_name = p_name.rsplit("/")
            f.write("Pin {}\n".format(inst_infos[inst_name]["alias"]+"/"+cell_infos["TA"][inst_infos[inst_name]["type"]]["pins"][pin_name]["alias"]))

    f.close()

def main(argv):
    design_name = ""
    top_name = ""
    try:
        opts, args = getopt.getopt(argv,"hn:t:",["name=","top="])
    except getopt.GetoptError:
        print('python lef2txt.py -n <design_name> -t <top>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('python lef2txt.py -n <design_name> -t <top>')
            sys.exit()
        elif opt in ("-n", "--name"):
            design_name = arg
        elif opt in ("-t", "--top"):
            top_name = arg

    parseLef("./pdks/Nangate15/back_end/lef/NanGate_15nm_OCL.macro.lef", "TA")
    parseLef("./pdks/Nangate45/Nangate45_stdcell.lef", "TB")
    d_x, d_y = parseDef("./datasets/"+design_name+"/Nangate15/"+top_name+".mapped.def")
    # pdb.set_trace()
    d_x, d_y = diesize(d_x, d_y)
    # parseDef("../../datasets/"+design_name+"/Nangate45/"+design_name+".mapped.def", "Nangate45")
    output("./place/input/"+design_name+".input", d_x, d_y)

if __name__ == "__main__":
   main(sys.argv[1:])