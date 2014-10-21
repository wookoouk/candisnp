#Installing CandiSNP
Navigate to a suitable web-root directory in your operating system and clone the repository:

`git clone candisnp`

###Dependancies
CandiSNP requires the following Perl modules to be installed:

* Statistics::R
* Tie::Handle::CSV
* Sort::Key::Natural
* Digest::MD5
* IPC::Open2
* Storable
* Number::Bytes::Human

CandiSNP also requires 'R' to be installed and we recommend installation of Java 1.8.
It also needs some R libraries to be installed:

* library(ggplot2) [0.9.3.1]
* library(gridExtra) [0.9.1]
* library(plyr) [1.8]
* library(ggthemes) [1.6.0]

CandiSNP is tested using R version 3.0.2 (2013-09-25).

When installing on a Linux-type environment, it might be necessary to edit /usr/bin/R to adjust the value of 'run_arch' to read as follows:
`run_arch=x86_64`

Also, it may be necessary to include a .htaccess file in both the candisnp base directory and the /public directory:
We suggest this file has the following text.

```
 AddType text/x-component htc
 AddType application/x-shockwave-flash swf
 AddType image/svg+xml svg
```
