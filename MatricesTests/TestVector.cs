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

        [Test]
        public void FindUnionWithVectorTest2()
        {
            Vector<double> vector1 = new Vector<double>(2);
            vector1[0] = 1;
            vector1[1] = 2;

            Vector<double> vector2 = new Vector<double>(2);
            vector2[0] = 2;
            vector2[1] = 3;

            Vector<double> resultVector = vector1.FindUnionWithVector(vector2);


            Vector<double> expectedVector = new Vector<double>(3);
            expectedVector[0] = 1;
            expectedVector[1] = 2;
            expectedVector[2] = 3;

            Assert.AreEqual(expectedVector.Length, resultVector.Length);

            Assert.AreEqual(expectedVector[0], resultVector[0]);
            Assert.AreEqual(expectedVector[1], resultVector[1]);
            Assert.AreEqual(expectedVector[2], resultVector[2]);
        }

        [Test]
        public void FindUnionWithVectorTest3()
        {
            Vector<double> vector1 = new Vector<double>(4);
            vector1[0] = 1;
            vector1[1] = 2;
            vector1[2] = 5;
            vector1[3] = 2;

            Vector<double> vector2 = new Vector<double>(1);
            vector2[0] = 1;

            Vector<double> resultVector = vector1.FindUnionWithVector(vector2);


            Vector<double> expectedVector = new Vector<double>(3);
            expectedVector[0] = 1;
            expectedVector[1] = 2;
            expectedVector[2] = 5;

            Assert.AreEqual(expectedVector.Length, resultVector.Length);

            Assert.AreEqual(expectedVector[0], resultVector[0]);
            Assert.AreEqual(expectedVector[1], resultVector[1]);
            Assert.AreEqual(expectedVector[2], resultVector[2]);
        }

        [Test]
        public void FindUnionWithVectorTest4()
        {
            Vector<double> vector1 = new Vector<double>(5);
            vector1[0] = 5;
            vector1[1] = 1;
            vector1[2] = 7;
            vector1[3] = 2;
            vector1[4] = 4;

            Vector<double> vector2 = new Vector<double>(1);
            vector2[0] = 8;

            Vector<double> resultVector = vector1.FindUnionWithVector(vector2);


            Vector<double> expectedVector = new Vector<double>(6);
            expectedVector[0] = 5;
            expectedVector[1] = 1;
            expectedVector[2] = 7;
            expectedVector[3] = 2;
            expectedVector[4] = 4;
            expectedVector[5] = 8;

            Assert.AreEqual(expectedVector.Length, resultVector.Length);

            Assert.AreEqual(expectedVector[0], resultVector[0]);
            Assert.AreEqual(expectedVector[1], resultVector[1]);
            Assert.AreEqual(expectedVector[2], resultVector[2]);
            Assert.AreEqual(expectedVector[3], resultVector[3]);
            Assert.AreEqual(expectedVector[4], resultVector[4]);
            Assert.AreEqual(expectedVector[5], resultVector[5]);
        }

        [Test]
        public void FindUnionWithVectorTest5()
        {
            Vector<double> vector1 = new Vector<double>(5);
            vector1[0] = 5;
            vector1[1] = 1;
            vector1[2] = 7;
            vector1[3] = 2;
            vector1[4] = 4;

            Vector<double> vector2 = new Vector<double>(0);

            Vector<double> resultVector = vector1.FindUnionWithVector(vector2);
            
            Vector<double> expectedVector = new Vector<double>(5);
            expectedVector[0] = 5;
            expectedVector[1] = 1;
            expectedVector[2] = 7;
            expectedVector[3] = 2;
            expectedVector[4] = 4;

            Assert.AreEqual(expectedVector.Length, resultVector.Length);

            Assert.AreEqual(expectedVector[0], resultVector[0]);
            Assert.AreEqual(expectedVector[1], resultVector[1]);
            Assert.AreEqual(expectedVector[2], resultVector[2]);
            Assert.AreEqual(expectedVector[3], resultVector[3]);
            Assert.AreEqual(expectedVector[4], resultVector[4]);
        }

        [Test]
        public void FindUnionWithVectorNullTest()
        {
            Vector<double> vector1 = new Vector<double>(5);
            vector1[0] = 5;
            vector1[1] = 1;
            vector1[2] = 7;
            vector1[3] = 2;
            vector1[4] = 4;

            Vector<double> vector2 = null;
            
            try
            {
                Vector<double> resultVector = vector1.FindUnionWithVector(vector2);
                Assert.Fail("Should have thrown an exception!");
            }
            catch (Exception e)
            {
                Assert.IsTrue(true);
            }

        }

        [Test]
        public void FindIntersectionWithVectorTest1()
        {
            Vector<double> vector1 = new Vector<double>(1);
            vector1[0] = 0;

            Vector<double> vector2 = new Vector<double>(1);
            vector2[0] = 0;

            Vector<double> resultVector = vector1.FindIntersectionWithVector(vector2);

            Vector<double> expectedVector = new Vector<double>(1);
            expectedVector[0] = 0;

            Assert.AreEqual(expectedVector.Length, resultVector.Length);

            Assert.AreEqual(expectedVector[0], resultVector[0]);
        }

        [Test]
        public void FindIntersectionWithVectorTest2()
        {
            Vector<double> vector1 = new Vector<double>(2);
            vector1[0] = 1;
            vector1[1] = 2;

            Vector<double> vector2 = new Vector<double>(2);
            vector2[0] = 2;
            vector2[1] = 3;

            Vector<double> resultVector = vector1.FindIntersectionWithVector(vector2);

            Vector<double> expectedVector = new Vector<double>(1);
            expectedVector[0] = 2;

            Assert.AreEqual(expectedVector.Length, resultVector.Length);

            Assert.AreEqual(expectedVector[0], resultVector[0]);
        }

        [Test]
        public void FindIntersectionWithVectorTest3()
        {
            Vector<double> vector1 = new Vector<double>(2);
            vector1[0] = 1;
            vector1[1] = 2;

            Vector<double> vector2 = new Vector<double>(2);
            vector2[0] = 3;
            vector2[1] = 4;

            Vector<double> resultVector = vector1.FindIntersectionWithVector(vector2);

            Vector<double> expectedVector = new Vector<double>(0);

            Assert.AreEqual(expectedVector.Length, resultVector.Length);
        }

        [Test]
        public void FindIntersectionWithVectorTestNull()
        {
            Vector<double> vector1 = new Vector<double>(2);
            vector1[0] = 1;
            vector1[1] = 2;

            Vector<double> vector2 = null;

            try
            {
                Vector<double> resultVector = vector1.FindIntersectionWithVector(vector2);
                Assert.Fail("Should have thrown an exception!");
            }
            catch (Exception e)
            {
                Assert.IsTrue(true);
            }
        }

        [Test]
        public void RemoveDuplicatesFindMultiplicityTest1()
        {
            Vector<double> vector = new Vector<double>(8);
            vector[0] = 0;
            vector[1] = 0;
            vector[2] = 0;
            vector[3] = 1;
            vector[4] = 2;
            vector[5] = 3;
            vector[6] = 3;
            vector[7] = 3;

            Vector<double>[] resultVector = vector.RemoveDuplicatesFindMultiplicity();

            Vector<double>[] expectedVector = new Vector<double>[2];

            expectedVector[0] = new Vector<double>(4);
            expectedVector[1] = new Vector<double>(4);

            expectedVector[0][0] = 0;
            expectedVector[0][1] = 1;
            expectedVector[0][2] = 2;
            expectedVector[0][3] = 3;

            expectedVector[1][0] = 0;
            expectedVector[1][1] = 2;
            expectedVector[1][2] = 2;
            expectedVector[1][3] = 2;

            Assert.AreEqual(expectedVector[0][0], resultVector[0][0]);
            Assert.AreEqual(expectedVector[0][1], resultVector[0][1]);
            Assert.AreEqual(expectedVector[0][2], resultVector[0][2]);
            Assert.AreEqual(expectedVector[0][3], resultVector[0][3]);

            Assert.AreEqual(expectedVector[1][0], resultVector[1][0]);
            Assert.AreEqual(expectedVector[1][1], resultVector[1][1]);
            Assert.AreEqual(expectedVector[1][2], resultVector[1][2]);
            Assert.AreEqual(expectedVector[1][3], resultVector[1][3]);
        }

        [Test]
        public void RemoveDuplicatesFindMultiplicityTest2()
        {
            Vector<double> vector = new Vector<double>(9);
            vector[0] = 0;
            vector[1] = 0;
            vector[2] = 0;
            vector[3] = 1;
            vector[4] = 1;
            vector[5] = 2;
            vector[6] = 3;
            vector[7] = 3;
            vector[8] = 3;

            Vector<double>[] resultVector = vector.RemoveDuplicatesFindMultiplicity();

            Vector<double>[] expectedVector = new Vector<double>[2];

            expectedVector[0] = new Vector<double>(4);
            expectedVector[1] = new Vector<double>(4);

            expectedVector[0][0] = 0;
            expectedVector[0][1] = 1;
            expectedVector[0][2] = 2;
            expectedVector[0][3] = 3;

            expectedVector[1][0] = 0;
            expectedVector[1][1] = 2;
            expectedVector[1][2] = 3;
            expectedVector[1][3] = 3;

            Assert.AreEqual(expectedVector[0][0], resultVector[0][0]);
            Assert.AreEqual(expectedVector[0][1], resultVector[0][1]);
            Assert.AreEqual(expectedVector[0][2], resultVector[0][2]);
            Assert.AreEqual(expectedVector[0][3], resultVector[0][3]);

            Assert.AreEqual(expectedVector[1][0], resultVector[1][0]);
            Assert.AreEqual(expectedVector[1][1], resultVector[1][1]);
            Assert.AreEqual(expectedVector[1][2], resultVector[1][2]);
            Assert.AreEqual(expectedVector[1][3], resultVector[1][3]);
        }

        [Test]
        public void RemoveDuplicatesFindMultiplicityTestNull()
        {
            Vector<double> vector = null;

            try
            {
                Vector<double>[] resultVector = vector.RemoveDuplicatesFindMultiplicity();
                Assert.Fail("Should have thrown an exception!");
            }
            catch (Exception e)
            {
                Assert.IsTrue(true);
            }
        }
    }
}
