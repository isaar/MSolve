# MSolve
Open source numerical solver for computational mechanics problems


#About MKL
1) Install the MKL bindings:
Compute.Net bindings for MKL, Cuda, etc
Project site: https://github.com/allisterb/Compute.NET
Needs .Net core 2.0
Usage (see site for current instructions): 
    1) Add the Compute.NET package feed to your NuGet package sources: https://www.myget.org/F/computedotnet/api/v2
    2) Install the bindings package into your project: Install-Package Compute.Bindings.IntelMKL.
    3) (Optional) Install the native library package into your project: Install-Package Compute.Winx64.IntelMKL.
2) After installing these packages and the dlls:
go to project -> Properties -> Build and set Platform target from "Any CPU" to "x64" 