oProject = oDesktop.GetActiveProject()
oDesign  = oProject.GetActiveDesign()

local_var_array = oDesign.GetVariables()

print local_var_array
with open("vars.csv", "w") as fd:
    for var in local_var_array:
        fd.write("%s\t%s\n" % (var, oDesign.GetVariableValue(var)))
