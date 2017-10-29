#!/bin/bash

HTML="plot.html"
FLAG="%U"
PROB="0.8"
TITLE="The Runtime of Lattice program in various conditions"

function build_html {
	cat << END_END
<html>
<head>
  <meta http-equiv='refresh' content='10' />
  <script type='text/javascript' src='https://www.google.com/jsapi'></script>
  <script type='text/javascript'>
    google.load('visualization', '1.1', {packages: ['corechart','line']});
    google.setOnLoadCallback(drawChart);

    function drawChart() {
      var options = {
        title: '$TITLE',
        width:  1000,
        height: 600,
		vAxis: {
          title: 'Time (Seconds)'
        },
        hAxis: {title: 'Size of Lattice'}
      };
      var chart = new google.visualization.LineChart(document.getElementById('linechart_material'));
      var data = new google.visualization.DataTable();
END_END
populatetest

cat << END_END
      chart.draw(data, options);
          }
  </script>
</head>
<body>
  <div id='linechart_material'></div>
</body>
</html>
END_END
}



function populatetest {

	echo "data.addColumn('string', 'Length of the lattice (Node x Node)');"
	echo "data.addColumn('number', 'Size of Site Lattice with Probability of 0.5 by Single Thread');"
	echo "data.addColumn('number', 'Size of Site Lattice with Probability of 0.5 by Multiple Threads');"
	echo "data.addColumn('number', 'Size of Bond Lattice with Probability of 0.5 by Single Thread');"
	echo "data.addColumn('number', 'Size of Bond Lattice with Probability of 0.5 by Multiple Thread');"
	for i in {256..1024..256}
	do
#time=$(TIMEFORMAT='%3R';time ./lab0 $i)

echo "data.addRows([['$i',$(/usr/bin/time -f "%e" ./lattice -size $i -p $PROB -s -l 2>&1),$(/usr/bin/time -f "%e" ./lattice -size $i -p $PROB -s 2>&1),
					  $(/usr/bin/time -f "%e" ./lattice -size $i -p $PROB -b -l 2>&1),$(/usr/bin/time -f "%e" ./lattice -size $i -p $PROB -b 2>&1)
					  ]]);"
	done



}



rm -f $HTML
build_html > $HTML
echo "output is in $HTML"

