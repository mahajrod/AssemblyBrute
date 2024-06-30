#!/usr/bin/env awk

{
if ( ($6 >0) && ($6<0.88)) print $0;
}
