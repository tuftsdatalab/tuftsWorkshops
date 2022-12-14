## Interactive Session

- Particularly good for debugging and working with software GUI. 

  **`$ srun [options] --pty [command]`**

- **Command** 

  - command to run an application, given the module is already loaded.
  - `bash` for a bash shell

- **Options**

  - Pseudo terminal `--pty`
  - Partition `-p` 
    - Default batch if not specified
    - You can start interactive sessions on any partition you have access to
  - Time `-t` or `--time=`
    - Default 15 minutes if not specified on non-"interactive" partition
  - Number of CPU cores `-n` 
    - Default 1 if not specified
  - Memory `--mem=`
    - Default 2GB if not specified
  - GPU `--gres=`
    - Default none
  - X Window `--x11=first`
    - Default none	

Starting an interactive session of bash shell on preempt partition with **2 CPU cores and 4GB of RAM, with X11 forwarding for 1 day, 2 hours, and 10 minutes** (use `exit` to end session and release resources).

```
[your_utln@login-prod-01 ~]$ srun -p batch --time=1-2:10:00 -n 2 --mem=4g --reservation=bioworkshop --pty bash
[your_utln@p1cmp025 ~]$
```

You will be placed on one of the reserved nodes for the workshop `p1cmp[025-026,028-029,043-045]`

The reservation will expire after the workshop. You will no longer have access to the reservation `bioworkshop`. 

In that case, you can simply omit the `--reservation=bioworkshop` option in the srun command

```
[your_utln@login-prod-01 ~]$ srun -p batch --time=1-2:10:00 -n 2 --mem=4g --pty bash

[your_utln@i2cmp003 ~]$ exit

```

