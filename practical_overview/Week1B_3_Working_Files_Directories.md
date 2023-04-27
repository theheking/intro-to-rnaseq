---
layout: page
title: Week 1A Intro to Shell - Navigating Files and Directories
---

Working with Files and Directories
==================================

> Overview
> --------
> 
> **Questions**
> 
> *   How can I view and search file contents?
>     
> *   How can I create, copy and delete files and directories?
>     
> *   How can I control who has permission to modify a file?
>     
> *   How can I repeat recently used commands?
>     
> 
> **Objectives**
> 
> *   View, search within, copy, move, and rename files. Create new directories.
>     
> *   Use wildcards (`*`) to perform operations on multiple files.
>     
> *   Make a file read only.
>     
> *   Use the `history` command to view and repeat recently used commands.
>     

Working with Files
------------------

### Our data set: FASTQ files

Now that we know how to navigate around our directory structure, let’s start working with our sequencing files. We did a sequencing experiment and have two results files, which are stored in our `data` directory.

### Wildcards

Navigate to your `data` directory:

    $ cd [scratch location]/data
    

We are interested in looking at the FASTQ files in this directory. We can list all files with the .fastq.gz extension using the command:

    $ ls *.fastq.gz
    
    
The `*` character is a special type of character called a wildcard, which can be used to represent any number of any type of character. Thus, `*.fastq` matches every file that ends with `.fastq`.

> - What is your output if run the command below?
>    $ ls *fastq.gz

> - Can you do another wildcard string that outputs the same as the command above `ls *fastq.gz`



`echo` is a built-in shell command that writes its arguments, like a line of text to standard output. The `echo` command can also be used with pattern matching characters, such as wildcard characters. Here we will use the `echo` command to see how the wildcard character is interpreted by the shell.

>     $ echo *.fastq.gz


Command History
---------------

If you want to repeat a command that you’ve run recently, you can access previous commands using the up arrow on your keyboard to go back to the most recent command. Likewise, the down arrow takes you forward in the command history.

A few more useful shortcuts:

*   Ctrl+C will cancel the command you are writing, and give you a fresh prompt.
*   Ctrl+R will do a reverse-search through your command history. This is very useful.
*   Ctrl+L or the `clear` command will clear your screen.

You can also review your recent commands with the `history` command, by entering:

    $ history
    

to see a numbered list of recent commands. You can reuse one of these commands directly by referring to the number of that command.

For example, if your history looked like this:

    259  ls *
    260  ls /srv/scratch/*.fastq.gz
    261  ls *R1*fastq.gz
    

then you could repeat command #260 by entering:

    $ !260
    

Type `!` (exclamation point) and then the number of the command from your history. You will be glad you learned this when you need to re-run very complicated commands. For more information on advanced usage of `history`, read section 9.3 of [Bash manual](https://www.gnu.org/software/bash/manual/html_node/index.html).

> Exercise
> --------
> 
> Find the line number in your history for the command that listed all the .fastq.gz files in your data directory. Rerun that command.
> 

Examining Files
---------------

We now know how to switch directories, run programs, and look at the contents of directories, but how do we look at the contents of files?

One way to examine a file is to print out all of the contents using the program `cat`.

Enter the following commands from within the `data` directory. First we want to unzip the fastq.gz to fastq

    $ gunzip SRR2589044_1.fastq.gz
    
    $ cat SRR2589044_1.fastq
    

This will print out all of the contents of the `SRR2589044_1.fastq` to the screen.

The print out is totally overwhelming. **Hint** Ctrl+C

To access just the final lines, we have to use the method called **piping**, where you pipe using `|` character. It is when the output of one command into another command. It is like doing step 1 and step 2 of a recipe all at once, rather than in sequence. 

The commands are `head` and `tail` and they let you look at the beginning and end of a file, respectively.

The `-n` option to either of these commands can be used to print the first or last `n` lines of a file.

    $ cat SRR2589044_1.fastq | head 
    
    $ cat SRR2589044_1.fastq | head -n 4



> Exercise
> --------
> 
> 1.  Print out the contents of the `[scratch loc]/data/SRR097977.fastq` file. What is the last line of the file? **Hint use the command `tail`**
> 2.  From your home directory, and without changing directories, use one short command to unzip all files in the `data` directory.
> 

`cat` is a terrific program, but when the file is really big, it can be annoying to use. The program, `less`, is useful for this case. `less` opens the file as read only, and lets you navigate through it. The navigation commands are identical to the `man` program.

Enter the following command:

    $ less SRR2589044_1.fastq
    

Some navigation commands in `less`:
| key   | action                 |
| ----- | ---------------------- |
| Space | to go forward          |
| b     | to go backward         |
| g     | to go to the beginning |
| G     | to go to the end       |
| q     | to quit                |

`less` also gives you a way of searching through files. Use the “/” key to begin a search. Enter the word you would like to search for and press `enter`. The screen will jump to the next location where that word is found.

**Shortcut:** If you hit “/” then “enter”, `less` will repeat the previous search. `less` searches from the current location and works its way forward. Scroll up a couple lines on your terminal to verify you are at the beginning of the file. Note, if you are at the end of the file and search for the sequence “CAA”, `less` will not find it. You either need to go to the beginning of the file (by typing `g`) and search again using `/` or you can use `?` to search backwards in the same way you used `/` previously.

For instance, let’s search forward for the sequence `TTTTT` in our file. You can see that we go right to that sequence, what it looks like, and where it is in the file. If you continue to type `/` and hit return, you will move forward to the next instance of this sequence motif. If you instead type `?` and hit return, you will search backwards and move up the file to previous examples of this motif.

> Exercise
> --------
> 
> What are the next three nucleotides (characters) after the first instance of the sequence quoted above?


Remember, the `man` program actually uses `less` internally and therefore uses the same commands, so you can search documentation using “/” as well!

There’s another way that we can look at files, and in this case, just look at part of them. This can be particularly useful if we just want to see the beginning or end of the file, or see how it’s formatted.

    
Creating, moving, copying, and removing
---------------------------------------

Now we can move around in the file structure, look at files, and search files. But what if we want to copy files or move them around or get rid of them? Most of the time, you can do these sorts of file manipulations without the command line, but there will be some cases (like when you’re working with a remote computer like we are for this lesson) where it will be impossible. You’ll also find that you may be working with hundreds of files and want to do similar manipulations to all of those files. In cases like this, it’s much faster to do these operations at the command line.

### Copying Files

When working with computational data, it’s important to keep a safe copy of that data that can’t be accidentally overwritten or deleted. For this lesson, our raw data is our FASTQ files. We don’t want to accidentally change the original files, so we’ll make a copy of them and change the file permissions so that we can read from, but not write to, the files.

First, let’s make a copy of one of our FASTQ files using the `cp` command.

Navigate to the `data` directory and enter:

    $ cp SRR2589044_1.fastq SRR2589044_1-copy.fastq
    $ ls -F
        

We now have two copies of the `SRR098026.fastq` file, one of them named `SRR098026-copy.fastq`. We’ll move this file to a new directory called `backup` where we’ll store our backup data files.
    

### Moving / Renaming

We can move files around using the command `mv`:

    $ mkdir backup
    $ mv SRR2589044_1-copy.fastq backup
    $ ls backup
    

The `mv` command is also how you rename files. Let’s rename this file to make it clear that this is a backup:

    $ cd backup
    $ mv SRR2589044_1-copy.fastq SRR2589044_1-backup.fastq
    $ ls
        

### File Permissions

We’ve now made a backup copy of our file, but just because we have two copies, it doesn’t make us safe. We can still accidentally delete or overwrite both copies. To make sure we can’t accidentally mess up this backup file, we’re going to change the permissions on the file so that we’re only allowed to read (i.e. view) the file, not write to it (i.e. make new changes).

View the current permissions on a file using the `-l` (long) flag for the `ls` command:

    $ ls -l
    
    -rw-r--r-- 1 helkin helkin 43332 Nov 15 23:02 SRR2589044_1-backup.fastq
    

The first part of the output for the `-l` flag gives you information about the file’s current permissions. There are ten slots in the permissions list. The first character in this list is related to file type, not permissions, so we’ll ignore it for now. The next three characters relate to the permissions that the file owner has, the next three relate to the permissions for group members, and the final three characters specify what other users outside of your group can do with the file. We’re going to concentrate on the three positions that deal with your permissions (as the file owner).

![Permissions breakdown](../assets/img/rwx_figure.svg)

Here the three positions that relate to the file owner are `rw-`. The `r` means that you have permission to read the file, the `w` indicates that you have permission to write to (i.e. make changes to) the file, and the third position is a `-`, indicating that you don’t have permission to carry out the ability encoded by that space (this is the space where `x` or executable ability is stored.

Our goal for now is to change permissions on this file so that you no longer have `w` or write permissions. We can do this using the `chmod` (change mode) command and subtracting (`-`) the write permission `-w`.

    $ chmod -w SRR2589044_1-backup.fastq
    $ ls -l 
   
    -r--r--r-- 1 helkin helkin 43332 Nov 15 23:02 SRR2589044_1-backup.fastq
    

### Removing

To prove to ourselves that you no longer have the ability to modify this file, try deleting it with the `rm` command:

    $ rm SRR2589044_1-backup.fastq
    

You’ll be asked if you want to override your file permissions:

    rm: remove write-protected regular file ‘SRR2589044_1-backup.fastq’? 
    

You should enter `n` for no. If you enter `n` (for no), the file will not be deleted. If you enter `y`, you will delete the file. This gives us an extra measure of security, as there is one more step between us and deleting our data files.

**Important**: The `rm` command permanently removes the file. Be careful with this command. It doesn’t just nicely put the files in the Trash. They’re really gone.

By default, `rm` will not delete directories. You can tell `rm` to delete a directory using the `-r` (recursive) option. Let’s delete the backup directory we just made.

Enter the following command:

    $ cd ..
    $ rm -r backup
    

This will delete not only the directory, but all files within the directory. If you have write-protected files in the directory, you will be asked whether you want to override your permission settings.

> Exercise
> --------
> 
> Starting in the `data` directory, do the following:
> 
> 1.  Make sure that you have deleted your backup directory and all files it contains.
> 2.  Create a backup of each of your FASTQ files using `cp`. (Note: You’ll need to do this individually for each of the two FASTQ files. We haven’t learned yet how to do this with a wildcard.)
> 3.  Use a wildcard to move all of your backup files to a new backup directory.
> 4.  Change the permissions on all of your backup files to be write-protected.
> 

> Key Points
> ----------
> 
> *   You can view file contents using `less`, `cat`, `head` or `tail`.
>     
> *   The commands `cp`, `mv`, and `mkdir` are useful for manipulating existing files and creating new directories.
>     
> *   You can view file permissions using `ls -l` and change permissions using `chmod`.
>     
> *   The `history` command and the up arrow on your keyboard can be used to repeat recently used commands.


Adapted from the Data Carpentry Intro to Command Line -shell genomics https://datacarpentry.org/shell-genomics/

Licensed under CC-BY 4.0 2018–2021 by The Carpentries  
Licensed under CC-BY 4.0 2016–2018 by [Data Carpentry](http://datacarpentry.org)
