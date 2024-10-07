import sys
arg_list = sys.argv
ids_name = arg_list[1]
r1_name = arg_list[2]
r2_name = arg_list[3]
with open(ids_name, "r") as ids:
    with open(r1_name, "a") as r1:
        with open(r2_name, "a") as r2:
            for line in ids:
                r1_line = line.strip("\n") + "/1\n"
                r2_line = line.strip("\n") + "/2\n"
                r1.write(r1_line)
                r2.write(r2_line)