using ISAAR.MSolve.Matrices;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatricesTests
{
    [TestClass]
    class TestVector
    {
        [TestMethod]
        public void SortAscendingTest1()
        {
            Vector<double> vector1 = new Vector<double>(3);
            vector1[0] = 7;
            vector1[0] = 8;
            vector1[0] = 9;
            vector1.SortAscending();

            Vector<double> expectedVector = new Vector<double>(3);
            expectedVector[0] = 7;
            expectedVector[1] = 8;
            expectedVector[2] = 9;

            Assert.AreEqual(vector1[0], expectedVector[0]);
            Assert.AreEqual(vector1[1], expectedVector[1]);
            Assert.AreEqual(vector1[2], expectedVector[2]);

        }

    }
}
