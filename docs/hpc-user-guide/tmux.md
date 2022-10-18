## Tmux
tmux is a terminal multiplexer. 
It lets you switch easily between several programs in one terminal, detach them (they keep running in the background) and reattach them to a different terminal.

**[Cheat Sheet]**(https://tmuxcheatsheet.com/)

- New Tmux Window `tmux new -s mysession`
- Detach it `CTRL+b d`
- List Sessions `tmux ls`
- Reattach Session `tmux a -t mysession`

### How to use tmux on Tufts HPC Cluster

#### Load tmux module
`[your_utln@login-prod-01 ~]$ module load tmux`
Make a note of the login nodename `login-prod-01` where your tmux session lives.

#### Start your tmux session
`[your_utln@login-prod-01 ~]$ tmux new -s mysession

#### Start your Interactive session inside the tmux session, and run your programs

#### Detach your tmux session OR lose connection...
`CTRL+b d`
#### Get your work session back
Log back in to the cluster or start a new terminal
If you are allocated on a **different** login node than where your tmux session lives. Simply do: 
`[your_utln@login-prod-03 ~]$ ssh login-prod-01`
`[your_utln@login-prod-01 ~]$ module load tmux`

#### Check tmux sessions
`[your_utln@login-prod-01 ~]$ tmux ls`

#### Pick your session to reattach it
`[your_utln@login-prod-01 ~]$ tmux a -t mysession`
