// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ true, true, true, false, true, false, false, false, true, false, true, false, true, false, true, true, false, false, true, false, true, false, true, true, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, true ];
var arrayMetadata    = [ [ "1", "OC4_F_B.host", "1.83458142926319", "F_host" ], [ "2", "OC4_F_B.sym", "0.360484399258279", "F_sym" ], [ "3", "OC5_H_B.host", "1.05367917496729", "H_host" ], [ "4", "OC5_H_B.sym", "0.631635682204692", "H_sym" ], [ "5", "OC9_C_B.host", "1.90352634440261", "C_host" ], [ "6", "OC9_C_B.sym", "0.551378026672332", "C_sym" ], [ "7", "OD4_C_B.host", "1.32158322850844", "C_host" ], [ "8", "OD4_C_B.sym", "0.453413512617458", "C_sym" ], [ "9", "OD5_F_B.host", "3.29259655864779", "F_host" ], [ "10", "OD5_F_B.sym", "0.814813852860716", "F_sym" ], [ "11", "OD6_H_B.host", "1.62750102120359", "H_host" ], [ "12", "OD6_H_B.sym", "0.343568107894447", "H_sym" ], [ "13", "OF7_C_B.host", "0.894510862617496", "C_host" ], [ "14", "OF7_C_B.sym", "0.507627048226183", "C_sym" ], [ "15", "OF8_F_B.host", "1.94427929010211", "F_host" ], [ "16", "OF8_F_B.sym", "0.322444399786401", "F_sym" ], [ "17", "OF9_H_B.host", "2.84704788648513", "H_host" ], [ "18", "OF9_H_B.sym", "0.530441230431767", "H_sym" ], [ "19", "OI1_C_B.host", "4.59517041510935", "C_host" ], [ "20", "OI1_C_B.sym", "0.958017869235603", "C_sym" ], [ "21", "OI2_F_B.host", "3.03304124514003", "F_host" ], [ "22", "OI2_F_B.sym", "0.758697666119144", "F_sym" ], [ "23", "OI3_H_B.host", "1.26812279284068", "H_host" ], [ "24", "OI3_H_B.sym", "0.252987376799962", "H_sym" ], [ "25", "OJ13_C_B.host", "3.49112120488666", "C_host" ], [ "26", "OJ13_C_B.sym", "0.72776621719813", "C_sym" ], [ "27", "OJ14_F_B.host", "5.22480066888067", "F_host" ], [ "28", "OJ14_F_B.sym", "0.693118441595341", "F_sym" ], [ "29", "OJ15_H_B.host", "1.1633516917651", "H_host" ], [ "30", "OJ15_H_B.sym", "0.623019325356849", "H_sym" ], [ "31", "OL6_C_B.host", "1.69561273738836", "C_host" ], [ "32", "OL6_C_B.sym", "0.665784319459367", "C_sym" ], [ "33", "OL7_F_B.host", "2.98782495681843", "F_host" ], [ "34", "OL7_F_B.sym", "0.664441913683274", "F_sym" ], [ "35", "OL8_H_B.host", "1.23739739145214", "H_host" ], [ "36", "OL8_H_B.sym", "0.355180011765369", "H_sym" ], [ "37", "OM1_C_B.host", "3.18441504972918", "C_host" ], [ "38", "OM1_C_B.sym", "0.711513871590235", "C_sym" ], [ "39", "OM2_F_B.host", "1.13321940449651", "F_host" ], [ "40", "OM2_F_B.sym", "0.418826337938269", "F_sym" ], [ "41", "OM3_H_B.host", "1.35942271150853", "H_host" ], [ "42", "OM3_H_B.sym", "1.15295159572438", "H_sym" ], [ "43", "OR7_C_B.host", "1.67623969611455", "C_host" ], [ "44", "OR7_C_B.sym", "0.45984374708935", "C_sym" ], [ "45", "OR8_F_B.host", "2.77028587925243", "F_host" ], [ "46", "OR8_F_B.sym", "0.59856688992799", "F_sym" ], [ "47", "OR9_H_B.host", "0.898139101214859", "H_host" ], [ "48", "OR9_H_B.sym", "0.343568107894447", "H_sym" ] ];
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
