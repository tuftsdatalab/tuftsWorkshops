# Overview

This tutorial provides step-by-step instructions on setting up SSH keys for secure and password-free access to Tufts HPC, enhancing both convenience and security.



# SSH Keys 

SSH keys provide a secure, password-free method to access remote systems, such as Tufts HPC. This tutorial will guide you through setting up an SSH key for simpler and safer connections.



# Steps to connect to Tufts HPC using SSH keys:

## 1. Generate SSH Key Pair

Start by creating a private and a public key on your local machine:

**Mac, Linux, and Windows Subsystem for Linux (WSL)**:

- Open the terminal. **Do not log in to Tufts server** 

- Run: `ssh-keygen`

- You will be prompted to enter a file name and a passphrase:

  - **Filename**: Hit Enter to use the default location (`~/.ssh/id_rsa`).

    - Your default location may be different. Here is an example:

      ```
      Generating public/private ed25519 key pair.
      Enter file in which to save the key (/Users/james/.ssh/id_ed25519): 
      ```

      If this is the case, then `.ssh/id_ed25519` is where your key file is. 

  - **Passphrase**: Set a passphrase for added security. **Important: Do not skip the passphrase. It protects your private key from unauthorized use if your local machine is compromised.** 

    ```
    Enter passphrase (empty for no passphrase): 
    ```

    Ex: You can try something like `tuftshpc` 

## 2. Copy Public Key to Cluster

Transfer your public key to the Tufts HPC to allow secure access:

- Run the following command:

  ```
  ssh-copy-id -i ~/.ssh/id_rsa.pub myusername@login.pax.tufts.edu
  ```

  **Remember: This is still on your local computer, do not log in to the server** 

- Replace `myusername` with your Tufts username.

- **Alternative**: If `ssh-copy-id` isn't available:

  ```
  cat ~/.ssh/id_rsa.pub | ssh myusername@login.pax.tufts.edu "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys"
  ```

- When prompted, enter your Tufts password and approve the Duo login notification.



## 3. Test SSH Connection

Verify that you can SSH from your local computer to the cluster without a Tufts password:

Execute: `ssh myusername@login.pax.tufts.edu`

You will see prompt something like this:

```
Enter passphrase for key '/Users/username/.ssh/id_ed25519': 
```

Or

```
Enter passphrase for key '/Users/username/.ssh/id_rsa': 
```

You should now be able to log in using just your passphrase, bypassing the need for your Tufts password and Duo approval. 

Ex: 

If your passphrase is set to `tuftshpc`, then just enter

```
Enter passphrase for key '/Users/username/.ssh/id_rsa': tuftshpc
```



## 4. Additional Steps for Custom Key Names or Locations

If your private key is not named `id_rsa` or located in the default directory, specify its location when connecting:

```
ssh -i path/to/my_private_key myusername@login.pax.tufts.edu
```

In our example:

`path/to/my_private_key` = `/Users/username/.ssh/id_rsa`

