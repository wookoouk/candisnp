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
When installing on a Linux-type environment, it might be necessary to edit /usr/bin/R to adjust the value of 'run_arch' to read as follows:
`run_arch=x86_64`
