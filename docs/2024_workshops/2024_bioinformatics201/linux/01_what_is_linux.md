# Introduction to Linux/Unix 

## What is Unix?

Unix is a powerful, multiuser, multitasking operating system that was originally developed in the 1960s and 1970s at **Bell Labs** by **Ken Thompson**, **Dennis Ritchie**, and others. It has since become one of the most influential operating systems in the history of computing, serving as the foundation for many modern operating systems, including Linux, macOS, and various BSDs. **Dennis Ritchie** is also the creator of C. 
<div style="display: flex; align-items: center;">
  <img src="https://www.redhat.com/sysadmin/sites/default/files/styles/full/public/2022-11/Ken_Thompson_Dennis_Ritchie_PDP-11_0.jpg?itok=D9Au5qnp" alt="Ken Thompson (sitting) and Dennis Ritchie at PDP-11" style="height:200px; margin-right: 10px;"/>

  <img src="https://upload.wikimedia.org/wikipedia/commons/0/0e/The_C_Programming_Language%2C_First_Edition_Cover.svg" alt="The C Programming Language - Wikipedia" style="height:200px;" />
</div>

## Key features
- **Multiuser**: Unix allows multiple users to access the system simultaneously, each with their own environment, files, and processes.
- **Multitasking**:  Unix is capable of running multiple processes at the same time. This means that it can handle several tasks concurrently, like running applications, performing background tasks, and processing commands.
- **Hierarchical File System**: Unix uses a hierarchical file system structure, where files are organized in directories (folders), starting from a root directory (`/`). 
- **Simple and Consistent Interface:** Unix provides a command-line interface (CLI) where users interact with the system by typing commands. **The Unix philosophy emphasizes small, simple, and modular commands that can be combined in scripts to perform complex tasks**.

## What is Linux?

Linux is a free, open-source, and Unix-like operating system kernel that was originally created by **Linus Torvalds** in 1991. Over time, Linux has grown into a full-fledged operating system used worldwide across various types of devices, from servers and desktop computers to smartphones and embedded systems.

<div style="display: flex; align-items: center;">
  <img src="https://fossbytes.com/wp-content/uploads/2016/08/linus-torvald-first-linux-email.png" alt="The First Linux Announcement" style="height:500px; margin-right: 10px;"/>

  <img src="https://media.licdn.com/dms/image/v2/D5612AQHNUVhyAZO6BA/article-inline_image-shrink_1500_2232/article-inline_image-shrink_1500_2232/0/1693106715409?e=1732752000&v=beta&t=qddDOg1IPPZp7IaXbaIDsa5lAJe-cOrz1t8dX9CiT4Q" alt="Linus Torvalds" style="height:500px;" />
</div>



<img src="https://fossbytes.com/wp-content/uploads/2016/08/linus-torvald-first-linux-email.png" alt="Linus Torvalds's Famous Email — The First Linux Announcement" width="1000">
<img src="https://media.licdn.com/dms/image/v2/D5612AQHNUVhyAZO6BA/article-inline_image-shrink_1500_2232/article-inline_image-shrink_1500_2232/0/1693106715409?e=1732752000&v=beta&t=qddDOg1IPPZp7IaXbaIDsa5lAJe-cOrz1t8dX9CiT4Q" alt="Linus Torvalds" width="500">

**Popular Linux Distributions:**

​	•	**Ubuntu:** A user-friendly distribution popular for desktop and server use, based on Debian.

​	•	**Fedora:** A cutting-edge distribution often used by developers and those who want the latest features.

​	•	**Debian:** Known for its stability and extensive software repositories, often used in server environments.

​	•	**CentOS/AlmaLinux/Rocky Linux:** Enterprise-grade distributions derived from Red Hat Enterprise Linux (RHEL).

​	•	**Arch Linux:** A rolling release distribution known for its simplicity and customization, aimed at advanced users.

​	•	**Kali Linux:** A distribution designed for penetration testing and security research.

<img src="../images/linux_distros.png" width="800">

## Shell

**Shell** is like a translator and bridge between you and the operating system's core, the **kernel**. It takes the commands you type and interprets them, telling the kernel what actions to perform. 

<img src="../images/shell.png" width="500">

[Previous: Overview ](00_overview.md)                                                                    

[Next: Files and File system](02_files.md)


