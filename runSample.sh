#!/usr/bin/env sh

./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.1 -o fb_indiana_0.1.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.2 -o fb_indiana_0.2.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.3 -o fb_indiana_0.3.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.4 -o fb_indiana_0.4.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.5 -o fb_indiana_0.5.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.6 -o fb_indiana_0.6.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.7 -o fb_indiana_0.7.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.8 -o fb_indiana_0.8.txt --zero-index --weak --normalize
./msparse -i ~/test_graphs/fb-indiana.edges -n 29732 --param-list 0.9 -o fb_indiana_0.9.txt --zero-index --weak --normalize