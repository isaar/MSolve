using ISAAR.MSolve.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;


namespace MatricesTests
{
    [TestFixture]
    class TestVector
    {
        [Test]
        public void SortAscendingTest1()
        {
            Vector<double> vector1 = new Vector<double>(3);
            vector1[0] = 7;
            vector1[1] = 8;
            vector1[2] = 9;
            vector1.SortAscending();

            Vector<double> expectedVector = new Vector<double>(3);
            expectedVector[0] = 7;
            expectedVector[1] = 8;
            expectedVector[2] = 9;

            Assert.AreEqual(expectedVector[0], vector1[0]);
            Assert.AreEqual(expectedVector[1], vector1[1]);
            Assert.AreEqual(expectedVector[2], vector1[2]);

        }

        [Test]
        public void SortAscendingTest2()
        {
            Vector<double> vector2 = new Vector<double>(3);
            vector2[0] = 8;
            vector2[1] = 9;
            vector2[2] = 7;
            vector2.SortAscending();

            Vector<double> expectedVector = new Vector<double>(3);
            expectedVector[0] = 7;
            expectedVector[1] = 8;
            expectedVector[2] = 9;

            Assert.AreEqual(expectedVector[0], vector2[0]);
            Assert.AreEqual(expectedVector[1], vector2[1]);
            Assert.AreEqual(expectedVector[2], vector2[2]);

        }

        [Test]
        public void SortAscendingTest3()
        {
            Vector<double> vector3 = new Vector<double>(4);
            vector3[0] = 8;
            vector3[1] = 9;
            vector3[2] = 7;
            vector3[3] = 9;
            vector3.SortAscending();

            Vector<double> expectedVector = new Vector<double>(4);
            expectedVector[0] = 7;
            expectedVector[1] = 8;
            expectedVector[2] = 9;
            expectedVector[3] = 9;

            Assert.AreEqual(expectedVector[0], vector3[0]);
            Assert.AreEqual(expectedVector[1], vector3[1]);
            Assert.AreEqual(expectedVector[2], vector3[2]);
            Assert.AreEqual(expectedVector[3], vector3[3]);
        }

        [Test]
        public void SortAscendingTest4()
        {
            Vector<double> vector4 = new Vector<double>(1);
            vector4[0] = 1;
            vector4.SortAscending();

            Vector<double> expectedVector = new Vector<double>(1);
            expectedVector[0] = 1;

            Assert.AreEqual(expectedVector[0], vector4[0]);
        }

        [Test]
        public void SortAscendingTest5()
        {
            Vector<double> vector5 = new Vector<double>(3);
            vector5[0] = -9;
            vector5[1] = -8;
            vector5[2] = -7;
            vector5.SortAscending();

            Vector<double> expectedVector = new Vector<double>(3);
            expectedVector[0] = -9;
            expectedVector[1] = -8;
            expectedVector[2] = -7;

            Assert.AreEqual(expectedVector[0], vector5[0]);
            Assert.AreEqual(expectedVector[1], vector5[1]);
            Assert.AreEqual(expectedVector[2], vector5[2]);
        }

        [Test]
        public void SortAscendingNullTest()
        {
            Vector<double> vector5 = null;

            try
            {
                vector5.SortAscending();
                Assert.Fail("Should have thrown an exception!");
            } catch (Exception e)
            {
                Assert.IsTrue(true);
            }
        }

        [Test]
        public void FindUnionWithVectorTest1()
        {
            Vector<double> vector1 = new Vector<double>(1);
            vector1[0] = 0;

            Vector<double> vector2 = new Vector<double>(1);
            vector2[0] = 0;

            Vector<double> resultVector=vector1.FindUnionWithVector(vector2);


            Vector<double> expectedVector = new Vector<double>(1);
            expectedVector[0] = 0;

            Assert.AreEqual(expectedVector[0], resultVector[0]);
        }
    }
}
