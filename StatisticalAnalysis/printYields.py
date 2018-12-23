import re
import math
from sys import argv, stdout, stderr, exit
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
import ROOT
ROOT.gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("-u", "--uncertainties", default=False, action="store_true", help="Report the uncertainties from the fit(s) too")
parser.add_option("-f", "--fit", default='pre', help="Report yields from before the fit")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

errors = False
if options.uncertainties: 
    errors = True

file = ROOT.TFile.Open(args[0]);
prefit = file.Get("norm_prefit")
fit_s = file.Get("norm_fit_s")
fit_b = file.Get("norm_fit_b")
if prefit == None: stderr.write("Missing fit_s in %s. Did you run FitDiagnostics in a recent-enough version of combine and with --saveNorm?\n" % file);
if fit_s  == None: raise RuntimeError, "Missing fit_s in %s. Did you run FitDiagnostics with --saveNorm?" % file;
if fit_b  == None: raise RuntimeError, "Missing fit_b in %s. Did you run FitDiagnostics with --saveNorm?" % file;

iter = fit_s.createIterator()
#Headline = "%-30s %-30s     pre-fit   signal+background Fit  bkg-only Fit"%("Channel","Process") if (prefit and errors) else "%-30s %-30s  signal+background Fit  bkg-only Fit"%("Channel","Process")
if prefit and errors :
 headrow  = ["Channel","Process","Pre-fit","S+B Fit","B-Only Fit"]
 headline = ("{:40} {:25} {:^25} {:^25} {:^25}").format(*headrow)
elif prefit: 
 headrow  = ["Channel","Process","Pre-fit","S+B Fit","B-Only Fit"]
 headline = ("{:40} {:25} {:>20} {:>20} {:>20}").format(*headrow)
else : 
 headrow = ["Channel","Process","S+B Fit","B-Only Fit"]
 headline = ("{:40} {:25} {:>20} {:>20}").format(*headrow)

line = "".join(["-" for i in range(len(headline))])
#print headline
#print line

procs = ['ggH','VH','ttH','qqZZ','ggZZ','zjets','qqH']
yields = {}
yields_err = {}
for ifs in ['4mu','4e','2e2mu','4l']:
  yields[ifs] = {}
  yields_err[ifs] = {}
  for iproc in procs:
    yields[ifs][iproc] = 0
    yields_err[ifs][iproc] = 0

while True:
    norm_s = iter.Next()
    if norm_s == None: break;
    norm_b = fit_b.find(norm_s.GetName())
    norm_p = prefit.find(norm_s.GetName()) if prefit else None
    # we have to replace any non-standard characters with "_" otherwise the matching will screw up 
    proc_chan_name = (norm_s.GetName()).replace(".","_").replace(":","_").replace(",","_")
    
    for ifs in ['4mu','4e','2e2mu','4l']:
      if(proc_chan_name.find(ifs) != -1):
	for iproc in procs:
	  if(proc_chan_name.find(iproc) != -1 or (iproc == 'VH' and (proc_chan_name.find('WH') != -1 or proc_chan_name.find('ZH') != -1))):
	    if(options.fit == 'pre'):
	      yields[ifs][iproc] += norm_p.getVal()
	      yields_err[ifs][iproc] += norm_p.getError()*norm_s.getError()
	    elif(options.fit == 'sb'):
	      yields[ifs][iproc] += norm_s.getVal()
	      yields_err[ifs][iproc] += norm_s.getError()*norm_s.getError()
	    elif(options.fit == 'b'):
	      yields[ifs][iproc] += norm_b.getVal()
	      yields_err[ifs][iproc] += norm_b.getError()*norm_s.getError()
	    else:
	      raise RuntimeError, "Non valid option. Expected pre/sb/b"
  
    m = re.match(r"(\w+)/(\w+)", proc_chan_name);
    if m == None: m = re.match(r"n_exp_(?:final_)?(?:bin)+(\.\w+)_proc_(\.\w+)", proc_chan_name);
    if m == None: raise RuntimeError, "Non-conforming object name %s" % norm_s.GetName()
    if norm_b == None: raise RuntimeError, "Missing normalization %s for background fit" % norm_s.GetName()
    if prefit and norm_p and errors:
        row = ["%-40s"%m.group(1), "%-25s"%m.group(2), "%10.3f +/- %-10.3f"%(norm_p.getVal(), norm_p.getError()), "%10.3f +/- %-10.3f"%(norm_s.getVal(), norm_s.getError()),"%10.3f +/- %-10.3f"%(norm_b.getVal(), norm_b.getError())]
	#print("{:<40} {:25} {:10} {:10} {:10}").format(*row)
        #print "%-30s %-30s % 7.3f +/- % 7.3f % 7.3f +/- % 7.3f  % 7.3f +/- % 7.3f" % 
    else:
        if norm_p and prefit:
            row = ["%-40s"%m.group(1), "%-25s"%m.group(2), "%10.3f"%(norm_p.getVal()), "%10.3f"%(norm_s.getVal()),"%10.3f"%(norm_b.getVal())]
	    #print("{:<40} {:25} {:>20} {:>20} {:>20}").format(*row)
            #print "%-30s %-30s %7.3f %7.3f %7.3f" % (m.group(1), m.group(2), norm_p.getVal(),  norm_s.getVal(),  norm_b.getVal())
        else:
            row = ["%-40s"%m.group(1), "%-25s"%m.group(2), "%10.3f"%(norm_s.getVal()),"%10.3f"%(norm_b.getVal())]
	    #print("{:<40} {:25} {:>20} {:>20}").format(*row)
            #print "%-30s %-30s %7.3f %7.3f" % (m.group(1), m.group(2), norm_s.getVal(), norm_b.getVal())

bkg = {'4mu':0,'4e':0,'2e2mu':0,'4l':0}
bkg_err = {'4mu':0,'4e':0,'2e2mu':0,'4l':0}
total4l = {'4mu':0,'4e':0,'2e2mu':0,'4l':0}
total4l_err = {'4mu':0,'4e':0,'2e2mu':0,'4l':0}
proc4l = {'ggH':0,'VH':0,'ttH':0,'qqZZ':0,'ggZZ':0,'zjets':0,'qqH':0}
proc4l_err = {'ggH':0,'VH':0,'ttH':0,'qqZZ':0,'ggZZ':0,'zjets':0,'qqH':0}
for ifs in ['4mu','4e','2e2mu']:
  print 'Channel', ifs
  for iproc in procs:
    print '%s = %.2f+/-%.2f' % (iproc,yields[ifs][iproc],math.sqrt(yields_err[ifs][iproc]))
    if(iproc != 'qqH'):
      bkg[ifs] += yields[ifs][iproc]
      bkg_err[ifs] += yields_err[ifs][iproc]
      bkg['4l'] += yields[ifs][iproc]
      bkg_err['4l'] += yields_err[ifs][iproc]
    total4l[ifs] += yields[ifs][iproc]
    total4l_err[ifs] += yields_err[ifs][iproc]
    total4l['4l'] += yields[ifs][iproc]
    total4l_err['4l'] += yields_err[ifs][iproc]
    proc4l[iproc] += yields[ifs][iproc]
    proc4l_err[iproc] += yields_err[ifs][iproc]
    
print '----------------- Totals ------------------------'
for iproc in procs:
  print '%s = %.2f+/-%.2f' % (iproc,proc4l[iproc],math.sqrt(proc4l_err[iproc]))

print ''
for ifs in ['4mu','4e','2e2mu','4l']:
  print 'Background %s = %.2f+/-%.2f' % (ifs,bkg[ifs],math.sqrt(bkg_err[ifs]))
  print 'S+B %s = %.2f+/-%.2f' % (ifs,total4l[ifs],math.sqrt(total4l_err[ifs]))
