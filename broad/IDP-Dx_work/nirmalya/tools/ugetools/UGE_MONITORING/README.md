# UGE_MONITORING 

## ugetop

ugetop provides a 'top'-like view of jobs that are queued on UGE.  

Run it like so:

    ./ugetop [refresh_seconds=30]


An example view is:


          USER           QUEUE         r        qw       Eqw        dr        dt
      nshoresh       [pending]         0      8197         0         0         0
       *bhaas*       [pending]         0      1542         0         0         0
        cychen       [pending]         0      6943         0         0         0
       bcleary       [pending]         0      2110         0         0         0
        jlevin       [pending]         0       342         0         0         0
        rkolde       [pending]         0       324         0         0         0
       *bhaas*           short       286         0         0         0         0
      nshoresh            long       262         0         0         0         0
        ljerby       [pending]         0       196         0         0         0
       gustavo       [pending]         0       106         0         0         0
        cychen           short       102         0         0         0         0
       bcleary            long       100         0         0         0         0
        lyping       [pending]         0        68         0         0         0
            li       [pending]         0        65         0         0         0
      mzekavat       [pending]         0        59         0         0         0
      rcdelros            long        50         0         0         0         0
       cfreije       [pending]         0         0         7         0         0
          phil       [pending]         0        41         0         0         0
      veerapen            long        40         0         0         0         0
       gustavo            long        38         0         0         0         0
        Totals   -------------      1322     20117         8         0         0



The top 20 largest jobs are shown in addition to any that you (the ${USER}) have running.  Your jobs will show up with '*' highlighting them.

Running jobs show up in the corresponding queue (eg. short, long).  Those not running are grouped as' [pending]'.


To stop 'ugetop', click ctrl-c



