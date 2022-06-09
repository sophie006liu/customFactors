from chimerax.core.commands import run
"""
This file contains code to produce one .defattr file
Then applies the attributes file to the display i.e.
"color byattribute kdHydrophobicity palette ^RdYlBu "
"""

def separateRosFromPDB(txt): 
  rosetta_start = txt.find("Nterm") 
  rosetta_lines = txt[rosetta_start-4:].splitlines() 

  ter_indices = [i+3 for i in range(rosetta_start-4) if txt.startswith("TER", i)]
  pdb_start = txt.find("ATOM")
  pdb_lines = txt[pdb_start:ter_indices[-1]].splitlines() 
  return rosetta_lines, pdb_lines
  
#get separate rosetta chains from a list of rosetta lines
#returns: a list of list of strings, one chain = list of strings
#one string is a residue
def getRosChainList(rosetta_lines): 
  rosetta_chain_list = []
  start = -1
  end = -1

  for i in range(len(rosetta_lines)):
    line = rosetta_lines[i]
    if "END_POSE_ENERGIES_TABLE" in line:
      break
    terms = line.split()
    residue = terms[0] 
  
    if "NtermProtein" in residue: #case we are starting a new amino acid
      start = i
    elif "CtermProtein" in residue:
      end = i
      rosetta_chain_list.append(rosetta_lines[start:end+1]) 
      start = i + 1
  
  return rosetta_chain_list
  
#iterate through and display the rosetta chains
def display_rosetta_chains(rosetta_chain_list):
  for i in rosetta_chain_list:
    for j in i:
      print(j[0:3], end = " ")
    print("")
  print("done printing rosetta chains")

#get separate PDB chains

def getPDBChainList(pdb_lines): 
  pdb_chain_list = []
  start = 0
  end = -1

  for i in range(len(pdb_lines)): #iterate through lines of PDB string list
    line = pdb_lines[i]
    if "TER" in line: #need to start new chain
      end = i 
      pdb_chain_list.append(pdb_lines[start:end]) #chain list does not include TER
      start = i + 1
      
  return pdb_chain_list

#iterate through and display the pdb chains
def printPDBChains(pdb_chain_list):
  for i in pdb_chain_list:
    for j in i: 
      terms = j.split()
      print(terms[3], end = " ")
    print("")
  print("done printing pdb ids")

#creates energyScore.defattr in the same folder as the source code (this file)
def createDeattrFile(filename):
  i = filename.rfind("/")
  rootFolderPath = filename[:i]
  attribute_file_path = rootFolderPath + "energyScore.defattr"
  f = open(attribute_file_path, "w") #if it doesn't exist it will create the file, if it aleady exists we make a new one 
  f.close()

  f = open(attribute_file_path, "a")
  header = """#\tEnergy Score scale from Rosetta
#\tMore positive means more INSERT TREND
#\tUse this file to assign the attribute in ChimeraX or Chimera.
# """ + filename + """
attribute: defattr
recipient: residues
"""
  f.write(header)
  f.close()
  return attribute_file_path

def writeDeattrFile(attribute_file_path, rosetta_chain_list, pdb_chain_list):
  f = open(attribute_file_path, "a+")
  for chain_index in range(len(rosetta_chain_list)):
    chain = rosetta_chain_list[chain_index]
    pdb_chain = pdb_chain_list[chain_index] #list of residues in string form for PDB 
    ros_chain = rosetta_chain_list[chain_index] #list of res in string form for PDB 
    pdb_res_row = pdb_chain[0].split() 
    pdb_res_num = pdb_res_row[5]
    chain_label = pdb_res_row[4]
    ros_res_stats = ros_chain[0].split() #get the first residue then find the _
    ros_res_and_num = ros_res_stats[0]
    
    i = ros_res_and_num.find("_")

    ros_res_num = ros_res_and_num[i+1:]
    shift = int(pdb_res_num) - int(ros_res_num)
  
    for residue_row in chain:
      stats = residue_row.split() #list of all the terms in the row  
  
      ros_res_and_num = stats[0]
      #if the res num is located later in start of chain
      if ros_res_and_num.count("_") > 1: 
        i = ros_res_and_num.rfind("_")
      else: 
        i = ros_res_and_num.find("_")
      ros_res_num = ros_res_and_num[i+1:]
      residue_number = int(ros_res_num) + shift

      total_energy = stats[-1]
      f.write("\t/"+chain_label+":" +str(residue_number) + "\t" + total_energy + "\n")
       
  return 

def makeAndSaveDeattrFile(filepath, rosetta_chain_list, pdb_chain_list):
  attribute_file_path = createDeattrFile(filepath)
  writeDeattrFile(attribute_file_path, rosetta_chain_list, pdb_chain_list)

def customColor(session, filepath):
    #open the file
    i = filepath.rfind("/")
    rootFolderPath = filepath[:i]
    attribute_file_path = rootFolderPath + "energyScore.defattr"

    with open(filepath) as f:
        txt = "".join(f.readlines())
  
    rosetta_lines, pdb_lines = separateRosFromPDB(txt) 
    rosetta_chain_list = getRosChainList(rosetta_lines)
    pdb_chain_list = getPDBChainList(pdb_lines)
     
    makeAndSaveDeattrFile(filepath, rosetta_chain_list, pdb_chain_list)
    
    open_attr_file = "open " + attribute_file_path

    color_cmd = "color byattribute defattr palette redblue"
    
    run(session, "open " + filepath)

    run(session, open_attr_file)

    run(session, color_cmd)
   
    f.close() 

def register_command(logger):
    from chimerax.core.commands import register, CmdDesc, StringArg

    desc = CmdDesc(required = [('filepath', StringArg)], 
                   synopsis='takes root directory and shows scores')
    register('customColor', desc, customColor, logger=logger)

register_command(session.logger)
 
print("bye") 