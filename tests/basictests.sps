* Encoding: UTF-8.
file handle data /name="c:/cc/misc2/extensions/r/stats_bayes_ttest/tests".
get file="data/sleep.sav".
dataset name sleep.

* one sample t test.
stats bayes ttest variables=extra testvalue=1.5 /options priorscale="medium 1 2".
stats bayes ttest variables=extra testvalue=1.5 /options posterior=yes.


* two sample test.
stats bayes ttest variables = extra group=group /options posterior=no.
stats bayes ttest variables = extra group=group/options posterior=yes.
stats bayes ttest variables = extra group=group/options posterior=yes priorscale="ultrawide".

STATS BAYES TTEST VARIABLES=extra GROUP=group  PAIRED=NO  
/OPTIONS POSTERIOR=YES ITERATIONS=3000 PRIORSCALE=".7071 1 1.414 2"
/SAVE WORKSPACE=RETAIN MODELFILE="c:\temp\modelfile.Rdata".

* paired test.  First restructure data.
SORT CASES BY ID group.
CASESTOVARS
  /ID=ID
  /INDEX=group
  /GROUPBY=VARIABLE.

stats bayes ttest variables = extra.1 extra.2 paired=yes.
stats bayes ttest variables = extra.1 extra.2 paired=yes /options posterior=yes.

