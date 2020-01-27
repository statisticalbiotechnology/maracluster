## Setting up a vagrant environment (Ubuntu 18.04)

Install vagrant and virtualbox:
```
sudo apt install vagrant virtualbox linux-headers-generic
```

The vagrant install might be outdated, so downloading the newest version from https://www.vagrantup.com/downloads.html might be better.

To compile for the platform of choice, use the `manager.sh` script in this folder:
```
usage: ./manager.sh 
                  [[-h]] [[-a]]
                  [[-b branch]]|[[-s source_directory]]
                  [[-r release_directory]]
                  -p ubuntu|centos|fedora|win32|win64|osx

If no branch and source_directory is provided, the source
code from which the sourcecode is checked out from will be used.
Make sure that Vagrant and VirtualBox are up to date.
  -h     prints this help page
  -a     keeps the vagrant box alive (i.e. do not call vagrant destroy)
```

To compile the Windows binaries you need a Windows virtual box with Visual Studio 2013/2015 and edit the concerning lines in `manager.sh` to point to this box. Also, you might need to install the winrm plugins for older versions of Vagrant (below 2.1):

```
vagrant plugin install winrm
vagrant plugin install winrm-fs
```

