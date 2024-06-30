#!/usr/bin/env awk

{
if ($6 >=0.99) print $0;
}
