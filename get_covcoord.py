import sys     
#Used from the command line
#Input should be basecov file from bbmap
#output is a csv file that contain contig name and start and stop of regions where reads have mapped (one coordinate per line)
# To run write in terminal: python get_covcoord.py input output_file                                                         
input_file = sys.argv[1]                                                      
output_file = sys.argv[2]                                              
scaffold = ""                                                            
found_start = False                                                      
start_coord = ""                                                         
coords = []             

prev_scaffold_name = ""
with open(input_file, "r") as basecov:                                        
    for line in basecov:                                                 
        if line[0] != "#":                                                                         
            split_line = line.strip("\n").split()                        
            scaffold_name = split_line[0].split()[0]   
            #If we are no longer looking at the same scaffold and previous coordinate have coverage, add last region of previous scaffold to list of regions with coverage
            if prev_scaffold_name != scaffold_name and found_start == True:
                coords.append([curr_scaff, str(start_coord), str(prev_coord)])
                found_start = False
            #If we have not found a start but there is coverage on the current coordinate, start a new region
            if int(split_line[-1]) != 0 and found_start == False:        
                found_start = True                                       
                start_coord = int(split_line[-2])
                curr_scaff = scaffold_name    
            #If we have found a start and there is no coverage on the current coordinate, append start coordinate and previous coordinate to coordinate list                        
            if int(split_line[-1]) == 0 and found_start == True:         
                found_start = False                                      
                stop_coord = prev_coord
                if scaffold_name == curr_scaff:                                    
                    coords.append([split_line[0], str(start_coord), str(stop_coord)])
                else:
                    coords.append([curr_scaff, str(start_coord), str(stop_coord)])
            #If we have found a start and there is coverage on the current coordinate, then update which coordinate was the previous one (for next iteration)                        
            if int(split_line[-1]) != 0 and found_start == True:
                prev_coord = int(split_line[-2]) 
            #Update which scaffold was the previous scaffold (for next iteration)
            prev_scaffold_name = scaffold_name
#When done, if applicable, add regions which were potentially at the end of the last scaffold
if int(split_line[-1]) != 0 and found_start == True:                     
    found_start = False                                                  
    stop_coord = int(split_line[-2])                                     
    coords.append([split_line[0], str(start_coord), str(stop_coord)])    
#Go through scaffolds and regions and write them to output                                                
with open(output_file, "a") as out_file:                                 
    for coord in coords:                                                                                                   
        out_file.write(coord[0] + "," + coord[1] + "," + coord[2] + "\n")
