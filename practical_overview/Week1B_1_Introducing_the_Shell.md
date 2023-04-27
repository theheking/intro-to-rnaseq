---
layout: page
title: Week 1A Intro to Shell - Introducting the Shell 
---

Introducing the Shell
=====================

> Overview
> --------
> 
> **Questions**
> 
> *   What is a command shell and why would I use one?
>     
> *   How can I move around on my computer?
>     
> *   How can I see what files and directories I have?
>     
> *   How can I specify the location of a file or directory on my computer?
>     
> 
> **Objectives**
> 
> *   Describe key reasons for learning shell.
>     
> *   Navigate your file system using the command line.
>     
> *   Access and read help files for `bash` programs and use help files to identify useful command options.
>     
> *   Demonstrate the use of tab completion, and explain its advantages.
>     

What is a shell and why should I care?
--------------------------------------

A _shell_ is a computer program that presents a command line interface which allows you to control your computer using commands entered with a keyboard instead of controlling graphical user interfaces (GUIs) with a mouse/keyboard/touchscreen combination.

There are many reasons to learn about the shell:

*   Bioinformatics team can only be used through command line, or have more features compared to the GUI.
*   Used for boring, repetitive tasks.
*   Ensures error-free output.
*   More reproducible output. 
*   Large amount of computational power and memory space is better. 


In this lesson you will learn how to use the command line interface to move around in your file system.

How to access the shell
-----------------------

On a Mac or Linux machine, you can access a shell through a program called “Terminal”, which is already available on your computer. The Terminal is a window into which we will type commands. If you’re using Windows, you’ll need to download a separate program to access the shell.

To save time, we are going to be working on a remote server where all the necessary data and software available. When we say a ‘remote sever’, we are talking about a computer that is not the one you are working on right now. This is Katana where the login instructions are Week 1A.

Clear Screen
-------------

This provides a lot of information about the remote server that you’re logging into. We’re not going to use most of this information for our workshop, so you can clear your screen using the `clear` command.

Type the word `clear` into the terminal and press the `Enter` key.

    $ clear
    

This will scroll your screen down to give you a fresh screen and will make it easier to read. You haven’t lost any of the information on your screen. If you scroll up, you can see everything that has been output to your screen up until this point.

> Tip
> ---
> 
> Hot-key combinations are shortcuts for performing common commands. The hot-key combination for clearing the console is `Ctrl+L`. Feel free to try it and see for yourself.

Navigating your file system
---------------------------

The part of the operating system that manages files and directories is called the **file system**. It organizes our data into files, which hold information, and directories (also called “folders”), which hold files or other directories.

Several commands are frequently used to create, inspect, rename, and delete files and directories.

Basic Commands - Where am I?
============================

The dollar sign is a **prompt**, which shows us that the shell is waiting for input; your shell may use a different character as a prompt and may add information before the prompt. When typing commands, either from these lessons or from other sources, do not type the prompt, only the commands that follow it.

Let’s find out where we are by running a command called `pwd` (which stands for “print working directory”). At any moment, our **current working directory** is our current default directory, i.e., the directory that the computer assumes we want to run commands in, unless we explicitly specify something else. Here, the computer’s response is `/home/[your_zID]`, which is the top level directory within our cloud system:

    $ pwd
  
Let’s look at how our file system is organized. We can see what files and subdirectories are in this directory by running `ls`, which stands for “listing”:

    $ ls
        

`ls` prints the names of the files and directories in the current directory in alphabetical order, arranged neatly into columns. 


> Basic Commands - Navigating to your directory 
> ----------------------------------------------

On the Katana HPC, you will have two locations: 

1. Home directory - which is the location where you are when you login
    - Small space, keep scripts or other small files here.
    
    
2. Scratch - where to keep large files 
    - Your scratch is /srv/scratch/[insert your zID here]
    - Large space, regularly cleaned of old files 

The command to change locations in our file system is `cd`, followed by a directory name to change our working directory. `cd` stands for “change directory”.

Let’s say we want to navigate to the directory we saw above. We can use the following command to get there:

    $ cd 
    
- Please navigate to your scratch space above.

We can make the `ls` output more comprehensible by using the **flag** `-F`, which tells `ls` to add a trailing `/` to the names of directories:

    $ ls -F
        
Anything with a “/” after it is a directory. Things with a “\*” after them are programs. If there are no decorations, it’s a file.

`ls` has lots of other options. To find out what they are, we can type:

    $ man ls

`man` (short for manual) displays detailed documentation (also referred as man page or man file) for `bash` commands. It is a powerful resource to explore `bash` commands, understand their usage and flags. Some manual files are very long. You can scroll through the file using your keyboard’s down arrow or use the Space key to go forward one page and the b key to go backwards one page. When you are done reading, hit q to quit.

> Challenge
> ---------
> 
> Use the `-l` option for the `ls` command to display more information for each item in the directory. What is one piece of additional information this long format gives you that you don’t see with the bare `ls` command? Which flag should you use to display human-readable format?

> > The additional information given includes the name of the owner of the file, when the file was last modified, and whether the current user has permission to read and write to the file.

No one can possibly learn all of these arguments, that’s what the manual page is for. You can (and should) refer to the manual page or other help files as needed.


### Shortcut: Tab Completion

Typing out file or directory names can waste a lot of time and it’s easy to make typing mistakes. Instead we can use tab complete as a shortcut. When you start typing out the name of a directory or file, then hit the Tab key, the shell will try to fill in the rest of the directory or file name.

Return to your home directory:

    $ cd
    

then enter:

    $ cd she<tab>
    

> Basic Commands - Make a Directory Downloading Trial Data
> ----------------------------------------------------------

Here we are using the -p option for mkdir. This option allows mkdir to create the new directory, even if one of the parent directories does not already exist. It also supresses errors if the directory already exists, without overwriting that directory.
 `wget` is short for “world wide web get”, and it’s basic function is to _download_ web pages or data at a web address.

It will take about 5 minutes to download the files.
**NB. Please make sure you are in your scratch directory **

    $   mkdir -p [yourscratch]/data/
    $   cd data

    $   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
    $   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz


These are two files with `.fastq.gz` extensions. FASTQ is a format for storing information about sequencing reads and their quality. We will be learning more about FASTQ files in a later lesson.

Using tab complete can be very helpful. However, it will only autocomplete a file or directory name if you’ve typed enough characters to provide a unique identifier for the file or directory you are trying to access.

For example, if we now try to list the files which names start with `SR` by using tab complete:

    $ ls SR<tab>
    


The shell auto-completes your command to `SRR2589044`, because all file names in the directory begin with this prefix. When you hit Tab again, the shell will list the possible choices.

    $ ls SRR09<tab><tab>
    
Tab completion can also fill in the names of programs, which can be useful if you remember the beginning of a program name.

    $ pw<tab><tab>

    pwck      pwconv    pwd       pwdx      pwunconv
    

Displays the name of every program that starts with `pw`.



Summary
-------

We now know how to move around our file system using the command line. This gives us an advantage over interacting with the file system through a GUI as it allows us to work on a remote server, carry out the same set of operations on a large number of files quickly, and opens up many opportunities for using bioinformatic software that is only available in command line versions.

In the next few episodes, we’ll be expanding on these skills and seeing how using the command line shell enables us to make our workflow more efficient and reproducible.

> Key Points
> ----------
> 
> *   The shell gives you the ability to work more efficiently by using keyboard commands rather than a GUI.
>     
> *   Useful commands for navigating your file system include: `ls`, `mkdir`, `wget`, `pwd`, and `cd`.
>     
> *   Most commands take options (flags) which begin with a `-`.
>     
> *   Tab completion can reduce errors from mistyping and make work more efficient in the shell.
>     

* * *
Adapted from the Data Carpentry Intro to Command Line -shell genomics  https://datacarpentry.org/shell-genomics/


