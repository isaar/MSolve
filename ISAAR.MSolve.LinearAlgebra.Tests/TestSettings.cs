using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests
{
    internal class TestSettings
    {
        public static TheoryData<LinearAlgebraProviderChoice> ProvidersToTest
            => new TheoryData<LinearAlgebraProviderChoice>
            {
                LinearAlgebraProviderChoice.Managed,
                LinearAlgebraProviderChoice.MKL
            };

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
