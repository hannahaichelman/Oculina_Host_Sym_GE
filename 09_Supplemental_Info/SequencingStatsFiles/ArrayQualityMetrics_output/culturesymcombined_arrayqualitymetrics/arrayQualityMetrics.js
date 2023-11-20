// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, true, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "OC4_F_B", "0.120269664590048", "freezing", "inhost", "freezing_inhost" ], [ "2", "OC5_H_B", "0.237437340955938", "hot", "inhost", "hot_inhost" ], [ "3", "OC9_C_B", "0.17252159674435", "control", "inhost", "control_inhost" ], [ "4", "OD4_C_B", "0.16090204653939", "control", "inhost", "control_inhost" ], [ "5", "OD5_F_B", "0.288817932795446", "freezing", "inhost", "freezing_inhost" ], [ "6", "OD6_H_B", "0.153154352199936", "hot", "inhost", "hot_inhost" ], [ "7", "OF7_C_B", "0.167300639202502", "control", "inhost", "control_inhost" ], [ "8", "OF8_F_B", "0.0822504014584964", "freezing", "inhost", "freezing_inhost" ], [ "9", "OF9_H_B", "0.193688036134687", "hot", "inhost", "hot_inhost" ], [ "10", "OI1_C_B", "0.388498680449581", "control", "inhost", "control_inhost" ], [ "11", "OI2_F_B", "0.298686715212052", "freezing", "inhost", "freezing_inhost" ], [ "12", "OI3_H_B", "0.0949454235819547", "hot", "inhost", "hot_inhost" ], [ "13", "OJ13_C_B", "0.308683356044934", "control", "inhost", "control_inhost" ], [ "14", "OJ14_F_B", "0.301912794302613", "freezing", "inhost", "freezing_inhost" ], [ "15", "OJ15_H_B", "0.294840844375586", "hot", "inhost", "hot_inhost" ], [ "16", "OL6_C_B", "0.274494989949209", "control", "inhost", "control_inhost" ], [ "17", "OL7_F_B", "0.266162908574188", "freezing", "inhost", "freezing_inhost" ], [ "18", "OL8_H_B", "0.140227281783112", "hot", "inhost", "hot_inhost" ], [ "19", "OM1_C_B", "0.310247062705503", "control", "inhost", "control_inhost" ], [ "20", "OM2_F_B", "0.119982373327085", "freezing", "inhost", "freezing_inhost" ], [ "21", "OM3_H_B", "0.491805110745392", "hot", "inhost", "hot_inhost" ], [ "22", "OR7_C_B", "0.175009342107538", "control", "inhost", "control_inhost" ], [ "23", "OR8_F_B", "0.190406708847874", "freezing", "inhost", "freezing_inhost" ], [ "24", "OR9_H_B", "0.099408852308739", "hot", "inhost", "hot_inhost" ], [ "25", "Control_1.1", "10.9738628384164", "control", "culture", "control_culture" ], [ "26", "Control_1", "10.4449388921122", "control", "culture", "control_culture" ], [ "27", "Control_2.1", "10.3050073875773", "control", "culture", "control_culture" ], [ "28", "Control_2", "6.21806405628357", "control", "culture", "control_culture" ], [ "29", "Control_3.1", "7.67837919745964", "control", "culture", "control_culture" ], [ "30", "Control_3", "4.62535285525981", "control", "culture", "control_culture" ], [ "31", "Control_4.1", "5.89044890465942", "control", "culture", "control_culture" ], [ "32", "Control_4", "4.71519981843034", "control", "culture", "control_culture" ], [ "33", "Cool_1", "3.46859298220203", "freezing", "culture", "freezing_culture" ], [ "34", "Cool_2", "5.44797011763534", "freezing", "culture", "freezing_culture" ], [ "35", "Cool_3", "8.45204818910703", "freezing", "culture", "freezing_culture" ], [ "36", "Cool_4", "9.35291193972726", "freezing", "culture", "freezing_culture" ], [ "37", "Heat_1.1", "8.71862437709667", "hot", "culture", "hot_culture" ], [ "38", "Heat_1", "9.40467216527592", "hot", "culture", "hot_culture" ], [ "39", "Heat_2.1", "8.01447321485116", "hot", "culture", "hot_culture" ], [ "40", "Heat_2", "9.52885931534103", "hot", "culture", "hot_culture" ], [ "41", "Heat_3.1", "9.72974786439434", "hot", "culture", "hot_culture" ], [ "42", "Heat_3", "4.89200596504868", "hot", "culture", "hot_culture" ], [ "43", "Heat_4.1", "9.93685138817417", "hot", "culture", "hot_culture" ], [ "44", "Heat_4", "5.60329702195336", "hot", "culture", "hot_culture" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
