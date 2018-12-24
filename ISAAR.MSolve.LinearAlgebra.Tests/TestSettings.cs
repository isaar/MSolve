using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests
{
    // Currently SuiteSparse dlls call MKL dll
    internal enum TestSuiteSparseAndMklLibs
    {
        None, MklOnly, MklAndSuiteSparse
    }

    internal class TestSettings
    {
        // Set the appropriate enums and flags here, in order to choose which native library tests will be run.
        internal static readonly TestSuiteSparseAndMklLibs librariesToTest = TestSuiteSparseAndMklLibs.None;

        internal const string MessageWhenSkippingMKL = "MKL is not set to be tested. See TestSettings.cs for more.";

        internal const string MessageWhenSkippingSuiteSparse
            = "SuiteSparse is not set to be tested. See TestSettings.cs for more.";

        public static TheoryData<LinearAlgebraProviderChoice> ProvidersToTest
        {
            get
            {
                var theoryData = new TheoryData<LinearAlgebraProviderChoice>();
                theoryData.Add(LinearAlgebraProviderChoice.Managed);
                if ((librariesToTest == TestSuiteSparseAndMklLibs.MklOnly) 
                    || (librariesToTest == TestSuiteSparseAndMklLibs.MklAndSuiteSparse))
                {
                    theoryData.Add(LinearAlgebraProviderChoice.MKL);
                }
                return theoryData;
            }
        }

        internal static bool TestMkl => (librariesToTest == TestSuiteSparseAndMklLibs.MklOnly) 
            || (librariesToTest == TestSuiteSparseAndMklLibs.MklAndSuiteSparse);

        internal static bool TestSuiteSparse => (librariesToTest == TestSuiteSparseAndMklLibs.MklAndSuiteSparse);

        internal static void RunMultiproviderTest(LinearAlgebraProviderChoice providers, Action test)
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
