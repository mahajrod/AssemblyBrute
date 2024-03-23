#!/usr/bin/env awk

{

SEQ_FIELD=10;
TAG_FIELD=12;
LEN=length($SEQ_FIELD);
NP=0;
RQ=0.0;

for (i=TAG_FIELD; i<=NF; i++) {
	if (substr($i,1,2) == "np") {
		split($i,NP_ARR,":");
		NP=NP_ARR[3];
		continue;
		};
	if (substr($i,1,2) == "rq") {
		split($i,RQ_ARR,":");
		RQ=RQ_ARR[3];
		continue;
		};
	};

printf "%s\t%i\t%.2f\t%.4f\t%i\t%.6f\n",$1,LEN,100*gsub(/C|G|c|g/,"",$SEQ_FIELD)/LEN,100*gsub(/N|n/,"",$SEQ_FIELD)/LEN,NP,RQ

}
