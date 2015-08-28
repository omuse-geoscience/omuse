import os, numpy, getopt, sys
import CdiObj


def usage():
  print("usage: " + sys.argv[0] + " [-o outputFile] [-e experiementTag] [-v] inputFile")

def htmlSkel(experiment):
  return ['<html> <head> <title>ECHAM Experiment: '+ experiment +'</title></head>' +
  """<body>
<script type="text/javascript" src="http://code.jquery.com/jquery-1.4.2.min.js"></script> 
<script type="text/javascript" src="https://code.zmaw.de/files/cdo/html/tableFold.js"></script>""",
  '''</body>
    </html>
  ''']

def tableSkel():
  return ["<table id=\"detail_table\" class=\"detail\"> <!--       <col style=\"width: 40px;\"> <col style=\"width: 80px;\"> <col style=\"width: 80px;\"> <col style=\"width: 80px;\">-->" +
          "<thead align=\"left\"><tr align=\"right\"><th>Code</th><th>North</th><th>South</th><th colspan=\"1000\">Global</th></tr></thead>",
          '</table>']

def tableRow(code,northVal,southVal,globalVal):
  return ["<tr  align=\"right\" title=\"Click to expand/collapse\" style=\"cursor: pointer;\" class=\"parent\" id=\"row%d\" >" % code +
          "<td valign=\"top\">%d</td><td>%3.3f</td><td>%3.3f</td><td>%3.3f</td>" % (code,northVal, southVal, globalVal) +
          '</tr>' + 
          "<tr style=\"display: none;\" class=\"child-row%d\">" % code + 
          '<td></td><td colspan="3">', 
          '</td></tr></table>']

def generateSubTable(headers, rows):
  html = []

  if len(headers) > 0:
    html.append("<tr>")
  for header in headers:
    html.append("<th>" + header + "</th>")
  html.append("</tr>")

  if len(rows) > 0:
    for row in rows:
      html.append("<tr>")
      for cell in row:
        html.append("<td>%3.3f</td>" % cell)
  html.append("</tr>")

  if html:
    html  = ["<table>"] + html + ["</table>"]

  return "".join(html)  

if __name__ == '__main__':
  try:
    opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
  except getopt.GetoptError, err:
    #print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

  output  = 'table.html'
  expName = '187'
  verbose = False
  for o, a in opts:
    if o == "-v":
      verbose = True
    elif o in ("-h", "--help"):
      usage()
      sys.exit()
    elif o in ("-o", "--output"):
      output = a
    elif o in ("-e", "--experiment"):
      expName = a
    else:
      assert False, "unhandled option"
  
  sourcefile = args[0]
  #TODO test if input file exists?
  if verbose == True:
    print 'Input:  ' + sourcefile
    print 'Output: ' + output
 
  cdoOperator = "zonavg"
  zonalfile   = cdoOperator + "-" + sourcefile
  cmd         = ' '.join(['cdo',cdoOperator, sourcefile, zonalfile])

  print cmd

  os.popen(cmd).read()

  cdi = CdiObj.Cdi(zonalfile)

  html = []
  htmlStart, htmlEnd = htmlSkel(expName)
  tableStart, tableEnd = tableSkel()

  html.append(htmlStart)
  html.append(tableStart)

  for code in cdi.codes:

      var = cdi.varByCode[code]
      var.getValues()
      vals = numpy.array(var.values)
      lats = numpy.array(var.grid.yvalues)
      northInd = numpy.where(lats >= 0)
      southInd = numpy.where(lats <  0)
      if verbose == True:
        print 'Global Mean:     ',vals.mean()
        print 'Global Max|Min:  ',vals.max(),'|',vals.min()
        print 'North  Mean:     ',vals[northInd].mean()
        print 'South  Mean:     ',vals[southInd].mean()
      
      codeTableStart, codeTableEnd = tableRow(code,vals[northInd].mean(), vals[southInd].mean(), vals.mean())
      
      html.append(codeTableStart)

      tdata = numpy.array([lats[northInd],vals[northInd],lats[southInd],vals[southInd]]).transpose().tolist()
      vtable  = generateSubTable(["latitude","values","latidute","values"],tdata)

      html.append(vtable)


  html.append(tableEnd)
  html.append(htmlEnd)

  f = open(output, 'w')
  f.write("".join(html))
  f.close
