using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Perhaps the ISubdomainDofOrderings should be accessed through IClusterDofOrdering. For now they are stored in Subdomain
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public interface IGlobalFreeDofOrdering
    {
        //TODO: Why do I need this ? Even if there were more than one subdomains, there would not 
        //      be a global vector (and obviously global matrix). Ideally all vectors should be on subdomain level, so that 
        //      they can be processed parallely (e.g. in a distributed environment).
        DofTable GlobalFreeDofs { get; }

        int NumGlobalFreeDofs { get; }

        IReadOnlyDictionary<ISubdomain_v2, ISubdomainFreeDofOrdering> SubdomainDofOrderings { get; }

        void AddVectorSubdomainToGlobal(ISubdomain_v2 subdomain, IVectorView subdomainVector, IVector globalVector);

        void AddVectorSubdomainToGlobalMeanValue(ISubdomain_v2 subdomain, IVectorView subdomainVector, IVector globalVector);

        void ExtractVectorSubdomainFromGlobal(ISubdomain_v2 subdomain, IVectorView globalVector, IVector subdomainVector);

        //TODO: the returned array should be readonly
        int[] MapFreeDofsSubdomainToGlobal(ISubdomain_v2 subdomain);
    }
}
