# Managing the environment

An environment is a collection of programs. More to come...

## Collection of tips

* How to invoke conda activate automatically on Windows powershell (JC, 2020-11-20)

    In modern windows, the default shell is powershell. When you open the terminal from VSCode, it starts the windows powershell but you can't run conda activate from there. To enable this, you have to run the following command as Administrator. Right click the Start Button, choose Windows Powershell (Admin) and run this command:
    
    > set-executionpolicy remotesigned

    After this, every time you open powershell, conda will automatically be activated on base (https://answers.microsoft.com/en-us/windows/forum/all/whats-wrong-with-my-windows-powershell/f05e72f2-a429-4ee0-81fb-910c8c8a1306?auth=1)