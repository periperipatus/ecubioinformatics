# Introduction to Unix

This tutorial is based on one by Dr. Nicholas Navin (MD Anderson Cancer Center, University of Texas). I’ve added a few minor things here and there based on what I think is most useful.

#Table of Contents

<!--ts-->
      *[How To](#how-to)
      *[Things to remember:](#things-to-remember)
      * [Introduction to Servers](#introduction-to-servers)
         * [Directories &amp; some basic unix terminology.](#directories--some-basic-unix-terminology)
            * [mkdir &amp; cd](#mkdir--cd)
         * [Jobs](#jobs)
         * [Doing things with files](#doing-things-with-files)
<!--te-->

### How To

If you are on Windows, you need to open your ```WSL``` (e.g. Ubuntu)... Or if you have selected  ```PuTTY``` (please don't) you'll need to login to the server first. 

If you are on Mac you'll just need to open your ```Terminal```

### Things to remember:
* Pay attention to the helpful hints here and there. They are… helpful
* Google is also helpful (seriously, use it when you are stumped, I do it all the time) 
* For any command, you can always ask for information about the usage and options using the inbuilt manual (man) page from the commandline. E.g.,: ```man grep``` or ```grep --help```
* Unix is CasE SensiTive
* It is very easy to overwrite files if you name them the same thing. (You will not be asked “are you sure”?)
* Keep your directories organized or things spiral out of control quickly.
* A good text editor is nice to have (on macs, TextWrangler is great, and free). Unix geeks like Unix text editors nano and vi that can be run in Unix. We are not going to cover those here.
* Even though I’ll show you how to transfer files using the command line, GUI ftp software are convenient (Cyberduck for macs, Filezilla, I think works on windows). (Again, a true Unix geek would frown on such things)
* you can use up and down arrows to scroll through previous commands 
* You can use ```tab``` to auto-complete names that are already in the system (This is useful for long filenames!)


<!-- toc -->

<!-- tocstop -->

## Introduction to Servers

Open your Mac Terminal or Linux WSL on Windows.

We are going to log into the university server using ```ssh```, which allows you to log into a remote computer. Please refer to your emails for the IP address and password. 

```bash
ssh ngsclass@<IP.ADDRESS>
```

You will be prompted for a password, type it in and click enter.


### Directories & some basic unix terminology.


To determine the directory you are in use ```pwd``` which stands for **Print Working Directory**

You are logged into the home directory. this can be accessed by the full name ```/home/ngsclass/``` or using ```~/``` or ```$HOME```

You can **list** the current directory contents using ```ls``` or ```ls .``` where the ```.``` indicates the current directory. There will be many files in the home directory.

There is also a root directory. To view what is in that directory you can use ```ls /``` 

The ```ls``` **command** also has a number of **arguments** that you can apply to it. For example ```ls -l```. Run it. Arguments in unix are typically preceded by a ```-```

**Question:** What does ```ls -l``` do? What do the different columns mean?

The output from a command is called **Standard Out** or **stdout**, and is usually output on your screen.

#### mkdir & cd 

Now, make a directory for you for the duration of the course using ```mkdir```. Ideally this will be your last name, and should not match any folders already existing in the home directory. I have made one for myself and called it ```Bolton/```

For example:

```bash
mkdir Xena
```

now move into the directory you just made (note: don't just copy Xena)

```bash
cd Xena
ls
```

You've made an empty directory! 

**Tip:** You can use ```cd``` to specify a full directory path that is not directly connected to your current directory (e.g. ```cd /usr/local/bin/``` is where most programs are installed).

You can also access the previous directory you were in using ```cd -```

You can jump back a directory using ```cd ../``` or two ```cd ../../``` and so on.

### Jobs 

Let’s run the sleep program on the server.  It will put the computer to sleep for a designated number of seconds

```bash
sleep 10
```

As you have noticed this runs the process in the “foreground” which can prevent you from continuing to issue commands.  Instead, use the ‘&’ sign at the end of the line to run the process in the “background”.


```bash
sleep 260 &
```

If you start a job in the **foreground**  (without using ```&```) and want to stop it you can kill it completely using ```ctrl+C``` or suspend it using ```ctrl+Z```.

Let's see what jobs your shell (session) has started. 

```bash
sleep 30
ctrl+z

jobs
```
	
You should be able to see two jobs - the ```[2] sleep 30``` and ```[1] sleep 260 &``` 
You can resume  ``` sleep 30``` and move it to the **background** by using ```bg %2```

We can see what processes are running using ```ps``` (**process status**)

You should still see ```sleep 260``` on the list of active processes. Use the Process ID (PID) in the first column to kill the process using the command ```kill <PID>```. Replacing everything including the brackets with the PID for the right proces.

Now, let's play with some genome data.

### Doing things with files

First, make sure you're in the directory with your name on it that we made at the beginning of the tutorial. 

you can make files using ```touch```

```bash
touch enzyme.txt
```

you can also make multiple files at once

```bash
 touch test test1 test12 test123 Test test.file test.file.txt
```

Note that you cannot make a filename with a space in it (it will think it is two separate files). Safe characters for use in filenames are ```a-zA-Z0-9._```

Now create a directory called ```bacteria```. 

**Move** the ```enzyme.txt``` file into this folder using ```mv```

```bash
mv enzyme.txt
cd bacteria
ls
```
you could also **copy** the file into the directory using ```cp```

Now you can **remove** the file using ```rm```.  Please be careful with this command **There is no recycle bin**. Removed files cannot be restored. 

```bash
rm enzyme.txt
```

**Question:** Look at ```man rm``` page for arguments you can pass to rm. Why might you want to run ```rm -i <FILE>```?

now, move back to the parent directory for ```bacteria``` (the one with your name on it). 

you can delete this now empty directory with

```bash
rmdir enzyme.txt
```

you may also use ```rm -r```for a full directory or ```rm -d``` removes only empty directories














