oProject = oDesktop.GetActiveProject()
oDesign = oProject.GetActiveDesign()

lengths = ["4.68mm", "5.50mm", "5.63mm", "5.63mm", "5.50mm", "4.68mm"]
widths  = ["4.56mm", "3.59mm", "3.39mm", "3.36mm", "3.39mm", "3.59mm", "4.56mm"]
wg_width = "7.112mm"
wg_height = "7.112mm/2"
wg_thickness = "2mm"
stub_length = "5mm"

def add_local_variable(name, value):
    oDesign.ChangeProperty(
	[
	    "NAME:AllTabs",
	    [
		"NAME:LocalVariableTab",
		[
		    "NAME:PropServers", 
		    "LocalVariables"
		],
		[
		    "NAME:NewProps",
		    [
			"NAME:%s" % name,
			"PropType:="		, "VariableProp",
			"UserDef:="		, True,
			"Value:="		, value
		    ]
		]
	    ]
	])

for i, length in enumerate(lengths):
    name = "length_%02i" % i
    add_local_variable(name, length)

for i, width in enumerate(widths):
    name = "width_%02i" % i
    add_local_variable(name, width)

add_local_variable("wg_width", wg_width)
add_local_variable("wg_height", wg_height)
add_local_variable("wg_thickness", wg_thickness)
add_local_variable("stub_length", stub_length)

    
oEditor = oDesign.SetActiveEditor("3D Modeler")
xpos = ["stub_length"]
for i, _ in enumerate(widths):
    x = '+'.join(xpos)
    oEditor.CreateBox(
	[
	    "NAME:BoxParameters",
	    "XPosition:="		, x,
	    "YPosition:="		, "width_%02i/2" % i,
	    "ZPosition:="		, "0mm",
	    "XSize:="		, "wg_thickness",
	    "YSize:="		, "(wg_width - width_%02i)/2" % i,
	    "ZSize:="		, "wg_height"
	], 
	[
	    "NAME:Attributes",
	    "Name:="		, "diaphragm_%02ip" % i,
	    "Flags:="		, "",
	    "Color:="		, "(132 132 193)",
	    "Transparency:="	, 0,
	    "PartCoordinateSystem:=", "Global",
	    "UDMId:="		, "",
	    "MaterialValue:="	, "\"pec\"",
	    "SolveInside:="		, False
	])
    oEditor.CreateBox(
	[
	    "NAME:BoxParameters",
	    "XPosition:="		, x,
	    "YPosition:="		, "-width_%02i/2" % i,
	    "ZPosition:="		, "0mm",
	    "XSize:="		, "wg_thickness",
	    "YSize:="		, "-(wg_width - width_%02i)/2" % i,
	    "ZSize:="		, "wg_height"
	], 
	[
	    "NAME:Attributes",
	    "Name:="		, "diaphragm_%02in" % i,
	    "Flags:="		, "",
	    "Color:="		, "(132 132 193)",
	    "Transparency:="	, 0,
	    "PartCoordinateSystem:=", "Global",
	    "UDMId:="		, "",
	    "MaterialValue:="	, "\"pec\"",
	    "SolveInside:="		, False
	])
    xpos.append("wg_thickness")
    if i < len(lengths):
        xpos.append("length_%02i" % i)

xpos.append("stub_length")
x = '+'.join(xpos)
print x
oEditor.CreateBox(
    [
	"NAME:BoxParameters",
	"XPosition:="		, 0,
	"YPosition:="		, "-wg_width/2",
	"ZPosition:="		, "0mm",
	"XSize:="		, x,
	"YSize:="		, "wg_width",
	"ZSize:="		, "wg_height"
    ], 
    [
	"NAME:Attributes",
	"Name:="		, "filter_body",
	"Flags:="		, "",
	"Color:="		, "(132 132 193)",
	"Transparency:="	, 0.9,
	"PartCoordinateSystem:=", "Global",
	"UDMId:="		, "",
	"MaterialValue:="	, "\"vacuum\"",
	"SolveInside:="		, True
    ])
