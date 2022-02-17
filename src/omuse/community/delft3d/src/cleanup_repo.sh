#!/bin/sh
find delft3dfm/ | grep '\.so$' | awk '{print "rm \""$0"\""}' | sh
find delft3dfm/ | grep '\.a$' | awk '{print "rm \""$0"\""}' | sh
find delft3dfm/ | grep '\.dll$' | awk '{print "rm \""$0"\""}' | sh
find delft3dfm/ | grep '\.o$' | awk '{print "rm \""$0"\""}' | sh
find delft3dfm/ | grep '\.exe$' | awk '{print "rm \""$0"\""}' | sh
find delft3dfm/ | grep '\.so\.' | awk '{print "rm \""$0"\""}' | sh
