# Starting with the Shell

We will spend most of our time learning about the basics of the shell by manipulating some experimental data that we download from the internet.

??? note "For Attendees Using Terminal Programs to Access the Cluster (instead of the Web Browser "OnDemand"

    If you are using a terminal on your home machine to connect to the tufts cluster, you will first need to log in by sending a simple command. **Ignore this if you are using the web browser login tool.**

    Replace "username01" with your tufts username.

    ```
    ssh username01@login.pax.tufts.edu
    ```

    Your username will have been created when your account was set up. If you do not have a cluster account, you can still follow this tutorial from your laptop or personal computer, except that the file structure will be different from what is described.

The login will ask you for your Tufts password.

??? note "Connection Issues?"

    If you are not on the Tufts network, you will need to set up the Tufts VPN (Virtual Private Network) before logging in:

    [VPN Instructions](https://it.tufts.edu/guides/vpn-virtual-private-network/anyconnect-desktop-application){:target="_blank" rel="noopener"}



!!! note "Best Practices for Logging In"

    If you are logged in to OnDemand, and on a machine called "login". If you are not on the login machine, type `exit` to get there.
    
    First, run the `tmux` commands (remember which login machine you are on, `login-prod-#`
    
    ```
    module load tmux
    tmux new -s OctoberWorkshop
    ```
    You will be able to recover this session if you are on the same login node and run `tmux a -t OctoberWorkshop`
    
    Second, run the `srun` command to go to a working machine, remember we have a reservation for this workshop on October 19 so if you are reading this at a different time, just drop the `--reservation=bioworkshop` parameter from the command.
    
    ```
    srun -p batch --time=1-2:10:00 -n 2 --mem=4g --reservation=bioworkshop --pty bash
    ```
The beginning of the line is called the 'command line prompt.'

```
[username01@login-prod-02 ~]$
```
It tells you who you are and what machine you have been assigned.

Once you run the `srun` command, the machine name should change

```
[username01@i2cmp003 ~]$
```

!!! question
    What machine are you on? Type your answer into the chat box.


!!! tip

    The `$` at the end of the line is where you start typing your commands. The `$` (on some Mac terminals it is a `%`) is not part of the command.
    The outputs from commands will not have that piece of information or `$` at the beginning of the line.


!!! tip

    The name of the computer you are on is important informatiom when troubleshooting the cluster. 
    `login` machines will reject large commands and output an error. 
    Make sure to switch machines with `srun` before running programs. 


## Using Basic Commands

Open up the shell through a terminal (OnDemand or on your laptop) and type the command::

```
whoami
```

and then hit ENTER 

(This is a good question for Mondays ....)

When you are on the Tufts cluster, this will return your username according to the cluster. This username is attached to you wherever you are in the cluster and creates a home where your files can be kept, regardless of which machine you are on in the cluster. [If you are on your laptop or personal computer, the answer to this may be different before you log in.]


## Running Commands

Let's try some simple commands.

Much like text shortcuts, shell commands often use abbreviations to get their point across.

For example, the command *pwd* is short for "print working directory." The word "print" here means it will output it into the visible screen.

Now type the command

```
pwd
```

You should see something similar to this:

```
/cluster/home/username01/
```

Try this command

```
ls
```

It may be empty for the moment, or it may not if this is not your first time using the shell. We will be creating content in the next part of this lesson.


!!! note "Takeaways"


    `pwd` and `ls` are examples of commands - programs you run at the shell
    prompt that do stuff. 

    * `pwd` stands for 'print working directory', while
    * `ls` stands for 'list files'. 

    It is similar to the abbreviations used in texting, it takes less time to get the point across (lol, tbh, imho, afaik, ftw -- you're saying them outloud in your head, right now, correct?)
