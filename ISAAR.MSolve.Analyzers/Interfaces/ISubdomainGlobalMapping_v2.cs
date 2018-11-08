//using ISAAR.MSolve.LinearAlgebra.Vectors;

////TODO: this should be hidden from analyzers and providers. Perhaps it will not be needed at all.
////TODO: Use IVectorView for the input vectors, at least for the global ones.
////TODO: Should I return the vectors, instead of void?
//namespace ISAAR.MSolve.Analyzers.Interfaces
//{
//    public interface ISubdomainGlobalMapping_v2
//    {
//        void SubdomainToGlobalVector(IVectorView subdomainVector, IVector globalVector);
//        void SubdomainToGlobalVectorMeanValue(IVectorView subdomainVector, IVector globalVector);
//        void SplitGlobalVectorToSubdomain(IVectorView globalVector, IVector subdomainVector);
//    }
//}
