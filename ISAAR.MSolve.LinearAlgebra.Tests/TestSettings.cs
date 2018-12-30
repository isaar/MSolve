using System;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests
{
    // Currently SuiteSparse dlls call MKL dll
    public enum TestSuiteSparseAndMklLibs
    {
        Neither, MklOnly, Both
    }

    public class TestSettings
    {
        // Set the appropriate enums and flags here, in order to choose which native library tests will be run.
        private static readonly TestSuiteSparseAndMklLibs librariesToTest = TestSuiteSparseAndMklLibs.Both;

        public const string MessageWhenSkippingMKL = "MKL is not set to be tested. See TestSettings.cs for more.";

        public const string MessageWhenSkippingSuiteSparse
            = "SuiteSparse is not set to be tested. See TestSettings.cs for more.";

        public static TheoryData<LinearAlgebraProviderChoice> ProvidersToTest
        {
            get
            {
                var theoryData = new TheoryData<LinearAlgebraProviderChoice>();
                theoryData.Add(LinearAlgebraProviderChoice.Managed);
                if ((librariesToTest == TestSuiteSparseAndMklLibs.MklOnly) 
                    || (librariesToTest == TestSuiteSparseAndMklLibs.Both))
                {
                    theoryData.Add(LinearAlgebraProviderChoice.MKL);
                }
                return theoryData;
            }
        }

        public static bool TestMkl => (librariesToTest == TestSuiteSparseAndMklLibs.MklOnly) 
            || (librariesToTest == TestSuiteSparseAndMklLibs.Both);

        public static bool TestSuiteSparse => (librariesToTest == TestSuiteSparseAndMklLibs.Both);

        public static void RunMultiproviderTest(LinearAlgebraProviderChoice providers, Action test)
        {
            LinearAlgebraProviderChoice defaultProviders = LibrarySettings.LinearAlgebraProviders; // Store it for later
            LibrarySettings.LinearAlgebraProviders = providers;

            try
            {
                test();
            }
            finally
            {
                LibrarySettings.LinearAlgebraProviders = defaultProviders; // Once finished, reset the default providers
            }
        }
    }
}
