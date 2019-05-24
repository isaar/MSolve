using System;
using System.Collections.Generic;
using System.Runtime.InteropServices.ComTypes;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Assemblers.Collocation;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class NurbsElement2DCollocation
	{
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = -1 ,Y= 0 , WeightFactor  = 1},
				new ControlPoint {ID = 1, X = -1.03750000000000 ,Y=  0,WeightFactor   = 1},
				new ControlPoint {ID = 2, X = -1.11250000000000, Y=   0, WeightFactor=   1},
				new ControlPoint {ID = 3, X = -1.22500000000000,Y=   0,WeightFactor   = 1},
				new ControlPoint {ID = 4, X = -1.37500000000000,Y=   0,WeightFactor   = 1},
				new ControlPoint {ID = 5, X = -1.56250000000000, Y=   0, WeightFactor   = 1},

				new ControlPoint {ID = 21, X = -1 ,Y= 0.00887131411164767,WeightFactor = 0.996338834764832},
				new ControlPoint {ID = 22, X = -1.03763434090205, Y=   0.00911916599990532, WeightFactor= 0.996430363895711},
				new ControlPoint {ID = 23, X = -1.11288260861762, Y=   0.00967068380654597,WeightFactor = 0.996610561872130},
				new ControlPoint {ID = 24, X = -1.22570594545548,Y=   0.0106654020184747,WeightFactor  = 0.996872277980737},
				new ControlPoint {ID = 25, X = -1.37605141302101 ,Y=  0.0123264255147491  ,WeightFactor = 0.997204071080174},
				new ControlPoint {ID = 26, X = -1.56385878147796,Y=   0.0149602450309121,WeightFactor  = 0.997590209601071},

				new ControlPoint {ID = 42,X= -0.999802541304798 ,Y= 0.0267260913477464, WeightFactor  = 0.989130915708095},
				new ControlPoint {ID = 43, X = -1.03769191515664 ,Y=  0.0274745622841582, WeightFactor  = 0.989402642815392},
				new ControlPoint {ID = 44, X = -1.11341015969418, Y=   0.0291367357095685, WeightFactor  = 0.989937605557885},
				new ControlPoint {ID = 45, X = -1.22684281257796, Y=   0.0321256995067823 ,WeightFactor = 0.990714575255314},
				new ControlPoint {ID = 46, X = -1.37783482437466,Y=   0.0371010766918943, WeightFactor  = 0.991699586019268},
				new ControlPoint {ID = 47, X = -1.56621120383911,Y=   0.0449674086989138,WeightFactor  = 0.992845934753180},

				new ControlPoint {ID = 63, X = -0.998902360107841  ,Y=0.0537346275140484,WeightFactor  = 0.978662271363786},
				new ControlPoint {ID = 64, X = -1.03713112918772, Y=   0.0552486885178388, WeightFactor  = 0.979195714579691},
				new ControlPoint {ID = 65, X = -1.11347091288754, Y=   0.0585990885122305, WeightFactor  = 0.980245930911005},
				new ControlPoint {ID = 66, X = -1.22770135613073,Y=   0.0645916607711512, WeightFactor  = 0.981771245106484},
				new ControlPoint {ID = 67, X = -1.37952691391903,Y=   0.0745107293408125,WeightFactor  = 0.983704976764141},
				new ControlPoint {ID = 68, X = -1.56861870501990,Y=   0.0901129591601654,WeightFactor  = 0.985955440331242},

				new ControlPoint {ID = 84, X = -0.996459496356937 ,Y= 0.0900905637327511,WeightFactor  = 0.965390547386301},
				new ControlPoint {ID = 85, X = -1.03505118290043, Y=   0.0926550612062497,WeightFactor  = 0.966255783701644},
				new ControlPoint {ID = 86, X = -1.11204726589193, Y=   0.0983007765271787, WeightFactor  = 0.967959217697474},
				new ControlPoint {ID = 87, X = -1.22710311280320,Y=    0.108320303774677 ,WeightFactor  = 0.970433252786657},
				new ControlPoint {ID = 88, X = -1.37976342932109,Y=   0.124768945837351 ,WeightFactor  = 0.973569734429773},
				new ControlPoint {ID = 89, X = -1.56953258980420,Y=   0.150448895513162  ,WeightFactor = 0.977219950135124},

				new ControlPoint {ID = 105, X = -0.991261829692497,Y=  0.135949159593355,WeightFactor   = 0.949945006550435},
				new ControlPoint {ID = 106, X = -1.03015194389203, Y=   0.139877140028209,WeightFactor   = 0.951196381386675},
				new ControlPoint {ID = 107, X = -1.10767121345326,Y=   0.148466439768619  ,WeightFactor = 0.953660025595520},
				new ControlPoint {ID = 108, X = -1.22335134114577,Y=   0.163553381883912 ,WeightFactor  = 0.957238175517892},
				new ControlPoint {ID = 109, X = -1.37658743092838,Y=   0.188050345547290  ,WeightFactor = 0.961774409299259},
				new ControlPoint {ID = 110, X = -1.56674302336333,Y=   0.225913844207604   ,WeightFactor = 0.967053646889642},
			};
		}

		private List<Knot> ElementKnot()
		{
			return new List<Knot>()
			{
				//new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
				//new Knot(){ID=1,Ksi=0.0,Heta=0.111111111,Zeta =0.0 },
				//new Knot(){ID=2,Ksi=0.111111111,Heta=0.0,Zeta =0.0 },
				//new Knot(){ID=3,Ksi=0.111111111,Heta=0.111111111,Zeta =0.0 }
			};
		}

		private LinearAlgebra.Vectors.Vector KnotValueVectorKsi()
		{
			return LinearAlgebra.Vectors.Vector.CreateFromArray(new double[46]
			{
				0, 0, 0, 0, 0, 0, 0.0312500000000000, 0.0625000000000000, 0.0937500000000000, 0.125000000000000,
				0.156250000000000, 0.187500000000000, 0.218750000000000, 0.250000000000000, .281250000000000,
				0.312500000000000, 0.343750000000000, 0.375000000000000, 0.406250000000000, .437500000000000,
				0.468750000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000,
				0.531250000000000, 0.562500000000000, 0.593750000000000, 0.625000000000000, 0.656250000000000,
				0.687500000000000, 0.718750000000000, 0.750000000000000, 0.781250000000000, 0.812500000000000,
				0.843750000000000, 0.875000000000000, 0.906250000000000, 0.937500000000000, 0.968750000000000, 1, 1, 1,
				1, 1, 1
			});

		}

		private LinearAlgebra.Vectors.Vector KnotValueVectorHeta()
		{
			return LinearAlgebra.Vectors.Vector.CreateFromArray(new double[27]
			{
				0, 0, 0, 0, 0, 0, 0.0625000000000000, 0.125000000000000, 0.187500000000000, 0.250000000000000,
				0.312500000000000, 0.375000000000000, 0.437500000000000, 0.500000000000000, 0.562500000000000,
				0.625000000000000, 0.687500000000000, 0.750000000000000, 0.812500000000000, 0.875000000000000,
				0.937500000000000, 1, 1, 1, 1, 1, 1
			});
		}

		private NURBSElement2DCollocation Element
		{
			get
			{
				//var element = new NURBSElement2D();
				NURBSElement2DCollocation element = new NURBSElement2DCollocation();
				var patch = new CollocationPatch();
				patch.Material = new ElasticMaterial2D(StressState2D.PlaneStrain)
				{
					YoungModulus = 100000,
					PoissonRatio = 0.3
				};


				foreach (var controlPoint in ElementControlPoints())
					element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
				foreach (var knot in ElementKnot())
					element.KnotsDictionary.Add(knot.ID, knot);
				patch.Thickness = 1;
				patch.DegreeKsi = 5;
				patch.DegreeHeta = 5;
				patch.NumberOfControlPointsHeta = 21;
				patch.KnotValueVectorKsi = KnotValueVectorKsi();
				patch.KnotValueVectorHeta = KnotValueVectorHeta();
				element.Patch = patch;
				element.CollocationPoint=new CollocationPoint2D(0,0.00625000000000000, 0.0125000000000000);
				return element;
			}
		}

        private NURBSElement2DCollocation BoundaryElement
        {
            get
            {
                //var element = new NURBSElement2D();
                NURBSElement2DCollocation element = new NURBSElement2DCollocation();
                var patch = new CollocationPatch();
                patch.Material = new ElasticMaterial2D(StressState2D.PlaneStrain)
                {
                    YoungModulus = 100000,
                    PoissonRatio = 0.3
                };


                foreach (var controlPoint in ElementControlPoints())
                    element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
                foreach (var knot in ElementKnot())
                    element.KnotsDictionary.Add(knot.ID, knot);
                patch.Thickness = 1;
                patch.DegreeKsi = 5;
                patch.DegreeHeta = 5;
                patch.NumberOfControlPointsHeta = 21;
                patch.KnotValueVectorKsi = KnotValueVectorKsi();
                patch.KnotValueVectorHeta = KnotValueVectorHeta();
                element.Patch = patch;
                element.CollocationPoint = new CollocationPoint2D(22, 0.00625000000000000, 0,true);
                return element;
            }
        }

        [Fact]
		private void TestCollocationNurbsValues()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var expectedValues = new double[36]
			{
				0.107756491046876, 0.172848409497250, 0.0446196272052037, 0.00352170471035475, 9.96040156364557e-05,
				8.76924569066375e-07, 0.172215582909445, 0.276270356282469, 0.0713301881900560, 0.00563137442349921,
				0.000159324600280913, 1.40325424039367e-06, 0.0441346527160369, 0.0708142973873237, 0.0182901074420221,
				0.00144472100564317, 4.09015227072497e-05, 3.60517704108796e-07, 0.00344655953090833,
				0.00553151875569626, 0.00142945498376100, 0.000112998516978540, 3.20222035571125e-06,
				2.82571933464716e-08, 9.61567751771517e-05, 0.000154379951879486, 3.99223781884040e-05,
				3.15901737750138e-06, 8.96349373609439e-08, 7.92114511665978e-10, 8.33030115505995e-07,
				1.33799464227073e-06, 3.46288694845191e-07, 2.74341649711324e-08, 7.79594672058808e-10,
				6.90131105601209e-12,
			};
			for (int i = 0; i < 36; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedValues[i], nurbs.NurbsValues[i, 0], 1e-9));
			}
		}


		[Fact]
		private void TestCollocationNurbsDerivativeValuesKsi()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var expectedDerivativeValues = new double[36]
			{
				-21.4902971135349, -34.4718321802253, -8.89866620951732, -0.702347300254383, -0.0198644171588461,
				-0.000174888485624470, 8.64765671172642, 13.8726772594869, 3.58178377490262, 0.282774601501794,
				0.00800034715608075, 7.04631993527534e-05, 11.2480889505601, 18.0476216978275, 4.66138833688865,
				0.368199348588330, 0.0104240984648629, 9.18809813721259e-05, 1.52909344026115, 2.45410211782634,
				0.634189027990329, 0.0501326872556841, 0.00142069042947834, 1.25365276876252e-05, 0.0600483608277138,
				0.0964077990130769, 0.0249308836131535, 0.00197275558582060, 5.59755779195311e-05, 4.94662782999498e-07,
				0.000666895671883862, 0.00107115315439964, 0.000277226990376325, 2.19628624949327e-05,
				6.24117067249048e-07, 5.52495568636568e-09,
			};
			for (int i = 0; i < 36; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i], nurbs.NurbsDerivativeValuesKsi[i, 0], 1e-9));
			}
		}

		[Fact]
		private void TestCollocationNurbsDerivativeValuesHeta()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var expectedDerivativeValues = new double[36]
			{
				-10.7764261575829, 4.28954547313487, 5.67289302078282, 0.780193880617294, 0.0310716423001000,
				0.000350763503953900, -17.2228001708165, 6.85614787894056, 9.06884597877821, 1.24756736468011,
				0.0497015803820970, 0.000561291577019773, -4.41378353511934, 1.75738469144463, 2.32538524761845,
				0.320061612345835, 0.0127592996624519, 0.000144204481880636, -0.344680806898117, 0.137274628716991,
				0.181739420717421, 0.0250335444667664, 0.000998938093238302, 1.13026735705029e-05, -0.00961637092281610,
				0.00383121734040580, 0.00507568966357562, 0.000699844600668407, 2.79617713550472e-05,
				3.16840092575466e-07, -8.33090186918566e-05, 3.32047536770791e-05, 4.40267846956413e-05,
				6.07772922226825e-06, 2.43195885572400e-07, 2.76047465571621e-09
			};
			for (int i = 0; i < 36; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i], nurbs.NurbsDerivativeValuesHeta[i, 0], 1e-9));
			}
		}


		[Fact]
		private void TestCollocationNurbsSecondDerivativeValuesKsi()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var expectedDerivativeValues = new double[36]
			{
				3423.75279424699, 5491.92182529329, 1417.70181859039, 111.895313455824, 3.16472375390534,
				0.0278625715677653, -4415.73302894191, -7083.77323668510, -1828.95799921366, -144.392543462572,
				-4.08519884146217, -0.0359804611907119, 552.395840640308, 886.322219103121, 228.921689739867,
				18.0823417720694, 0.511929507290087, 0.00451229290299885, 411.533530845430, 660.486326742618,
				170.682865441969, 13.4924925147660, 0.382358417927473, 0.00337402631388489, 27.6278335032534,
				44.3565583278856, 11.4705262901248, 0.907651135124908, 0.0257539743915116, 0.000227590908737324,
				0.427265525131462, 0.686264485294914, 0.177613291879907, 0.0140711274235435, 0.000399858205299781,
				3.53971871791393e-06,
			};
			for (int i = 0; i < 36; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i], nurbs.NurbsSecondDerivativeValueKsi[i, 0], 1e-9));
			}
		}

		[Fact]
		private void TestCollocationNurbsSecondDerivativeValuesKsiHeta()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var expectedDerivativeValues = new double[36]
			{
				2149.06069562792  ,
				-855.679659372179 ,
				-1131.41851180411 ,
				-155.601159440793 ,
				-6.19685336718003 ,
				-0.0699551420991173,
				-865.026185216468 ,
				343.957661160148 ,
				455.302192175189  ,
				62.6390378398084 ,
				2.49553859812347  ,
				0.0281831570096359,
				-1124.94078926736 ,
				447.802849169456 ,
				592.623022738136 ,
				81.5687424249061 ,
				3.25176820603158 ,
				0.0367513151596485, 
				-152.924336283433 ,
				60.8966118912998 ,
				80.6284890246988  ,
				11.1062003057394 ,
				0.443183114335965 ,
				0.00501448815852712, 
				-6.00537974270124 ,
				2.39235602803926 ,
				3.16964069157811  ,
				0.437038071570744  ,
				0.0174615737210826 ,
				0.000197860634670657, 
				-0.0666953348003726,
				0.0265810606788613 ,
				0.0352459529335453 ,
				0.00486559242236938 ,
				0.000194693484964216 ,
				2.20993449225890e-06,
			};
			for (int i = 0; i < 36; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i], nurbs.NurbsSecondDerivativeValueKsiHeta[i, 0], 1e-9));
			}
		}


		[Fact]
		private void TestCollocationNurbsSecondDerivativeValuesHeta()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var expectedDerivativeValues = new double[36]
			{
				862.208137050026 ,
				-1110.45756959740 , 
				136.328997697885  ,
				104.674079793138  ,
				7.13654935955457  ,
				0.112241292397297 ,
				1377.97431476070  ,
				-1774.88765607716 ,
				217.939361457008  ,
				167.378864564433  ,
				11.4154822657559  ,
				0.179608457853405 ,
				353.141201313712  ,
				-454.943570485826 ,
				55.8828518196742  ,
				42.9408068708527  ,
				2.93056192379451  ,
				0.0461441889857859 ,
				27.5774725355872   ,
				-35.5370170399551  ,
				4.36749873086471  ,
				3.35860520842093  ,
				0.229436569225423 ,
				0.00361675794316606,
				0.769393594620336  ,
				-0.991807715542706 ,
				0.121977213729529 ,
				0.0938940837567214 ,
				0.00642227274402007,
				0.000101386093687283, 
				0.00666544852211772 ,
				-0.00859589210517978 ,
				0.00105803642117494,
				0.000815413616253236, 
				5.58573449277392e-05 ,
				8.83328052933049e-07
			};
			for (int i = 0; i < 36; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i], nurbs.NurbsSecondDerivativeValueHeta[i, 0], 1e-9));
			}
		}


		[Fact]
		private void TestCollocation2DJacobian()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var jacobian=Element.CalculateJacobianMatrix(Element, nurbs);

			var expectedJacobian = new double[2, 2]
			{
				{ -0.00796405057272683  ,  1.46104086310680},
				{ -3.01022037002596 ,  0.0216716120925286}
			};
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedJacobian[i, j], jacobian[i, j], 1e-9));
				}
			}
		}

		[Fact]
		private void TestCollocation2DHessian()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var hessianMatrix = Element.CalculateHessian(Element, nurbs, 0);
			var expectedHessian = new double[3, 2]
			{
				{2.16476793547675 ,  0.852779461353878},
				{0.0426696295135573, 0.141574351233695 },
				{-1.60009113818852 , 3.47799425599668}
			};

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedHessian[i, j], hessianMatrix[i, j], 1e-9));
				}
			}
		}

		[Fact]
		private void TestCollocation2DNaturalDerivatives()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var jacobianMatrix = Element.CalculateJacobianMatrix(Element, nurbs);

			var inverseJacobian = jacobianMatrix.Invert();
			var dR = Element.CalculateNaturalDerivatives(nurbs, inverseJacobian);

			var expectedDR = new double[2, 36]
			{
				{3.47418786888943,-1.59491791717771, -1.92846834837001 , -0.262652804419645 ,-0.0104203407760361 ,-0.000117390572122759, 5.76427956322017 ,-2.20935170808089 ,-2.99515321425408 ,-0.413066690606795 ,-0.0164721684944232 ,-0.000186122049081198, 1.52175109688758  ,-0.494894952798661 ,-0.749556897897653 ,-0.104514760021170 ,-0.00418745878312805, -4.74540731078801e-05, 0.122042978389101 ,-0.0335114676486208,-0.0572513769355071,-0.00806946895051393,-0.000324861058261046, -3.69313678800837e-06 ,0.00349060169214760,-0.000797714074348734, -0.00156336549817325  ,-0.000222777392010279 ,-9.01347652676035e-06 ,-1.02821345353076e-07 ,3.09627619179099e-05 ,-5.75274283568431e-06 ,-1.32602399879104e-05 ,-1.91088332929543e-06  ,-7.77177456175474e-08 ,-8.89844530737662e-10},
				{ -14.6899577196008,-23.6027171162401,-6.10114669210320 ,-0.482148787388404 ,-0.0136528729506954 ,-0.000120341185874666,5.95025364807184 ,9.48302146819428 ,2.43520240471957 ,0.191291650045700 ,0.00538599718314209,4.72135418543246e-05,7.70697694882454 ,12.3498806809734  ,3.18637140507299  ,0.251441966498285 ,0.00711188139478022,6.26286759286379e-05,1.04724339705292  ,1.67951170481660  ,0.433754517845223 ,0.0342690083906348,0.000970611606690962,8.56041447935826e-06,0.0411187405315767 ,0.0659813482374857 ,0.0170552607531285 ,0.00134902549625684,3.82629915071373e-05, 3.38008279628147e-07 ,0.000456621219660495 ,0.000733112513353756 ,0.000189673945576838 ,1.50219235325079e-05 ,4.26749268235201e-07 ,3.77666981044437e-09}
			};

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 36; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedDR[i, j], dR[i, j], 1e-9));
				}
			}
		}

		[Fact]
		private void TestCollocation2DSquareDerivatives()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var jacobianMatrix = Element.CalculateJacobianMatrix(Element, nurbs);

			var sq = Element.CalculateSquareDerivatives(jacobianMatrix);

			var expectedSquareDerivatives = new double[3, 3]
			{
				{ 6.34261015249505e-05  ,  2.13464040366785 ,   -0.0232716066452060},
				{ 9.06142667611920  ,  0.000469658770689033 ,   -0.130472656344461},
				{ 0.0239735472619392 , 0.0316631108365837,  -4.39822756137908}
			};

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedSquareDerivatives[i, j], sq[i, j], 1e-9));
				}
			}
		}

		[Fact]
		private void TestCollocation2DSquareDerivativesInverse()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var jacobianMatrix = Element.CalculateJacobianMatrix(Element, nurbs);

			var sq = Element.CalculateSquareDerivatives(jacobianMatrix);

			var sqInv = sq.Invert();
			var expectedSquareDerivatives = new double[3, 3]
			{
				{ 2.42826017963574e-05  ,  0.110366560012572 ,  -0.00327413330725496},
				{ 0.468499747934614 ,  3.27929736000920e-06 ,   -0.00247899171968521},
				{ 0.00337289086997935, 0.000601601835837063  ,  -0.227400008576538}
			};

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedSquareDerivatives[i, j], sqInv[i, j], 1e-9));
				}
			}
		}



		[Fact]
		private void TestCollocation2DddR3()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var jacobianMatrix = Element.CalculateJacobianMatrix(Element, nurbs);
			var squareDerivatives = Element.CalculateSquareDerivatives(jacobianMatrix);

			var inverseJacobian = jacobianMatrix.Invert();
			var dR = Element.CalculateNaturalDerivatives(nurbs, inverseJacobian);

			var hessianMatrix = Element.CalculateHessian(Element, nurbs, 0);

			var ddR2 = new double[3, nurbs.NurbsSecondDerivativeValueKsi.NumRows];
			for (int i = 0; i < ddR2.GetLength(1); i++)
			{
				ddR2[0, i] = hessianMatrix[0, 0] * dR[0, i] + hessianMatrix[0, 1] * dR[1, i];
				ddR2[1, i] = hessianMatrix[1, 0] * dR[0, i] + hessianMatrix[1, 1] * dR[1, i];
				ddR2[2, i] = hessianMatrix[2, 0] * dR[0, i] + hessianMatrix[2, 1] * dR[1, i];
			}

			var ddR3 = new double[3, nurbs.NurbsSecondDerivativeValueKsi.NumRows];
			for (int i = 0; i < ddR2.GetLength(1); i++)
			{
				ddR3[0, i] = nurbs.NurbsSecondDerivativeValueKsi[i, 0] - ddR2[0, i];
				ddR3[1, i] = nurbs.NurbsSecondDerivativeValueHeta[i, 0] - ddR2[1, i];
				ddR3[2, i] = nurbs.NurbsSecondDerivativeValueKsiHeta[i, 0] - ddR2[2, i];
			}

			var expectedDDR3 = new double[3, 36]
			{
				{
					3428.75927797803,
					5515.50236484899 ,
					1427.07943762525  ,
					112.875062408196 ,
					3.19892426313487 ,
					0.0282193194058928,
					-4433.28561061242 ,
					-7087.07742888891,
					-1824.55087816859 ,
					-143.661479525820 ,
					-4.05413348705436,
					-0.0356178128855668,
					542.529251008904 ,
					876.861827233464  ,
					227.827034387759  ,
					18.0941674285951 ,
					0.514929517409880 ,
					0.00456161151034683,
					410.376268459018 ,
					659.126618206284 ,
					170.436904442933 ,
					13.4807371358899 ,
					0.382233949086744 ,
					0.00337472095233623,
					27.5852119432326 ,
					44.3020176551260 ,
					11.4593662375484  ,
					0.906982975243831 ,
					0.0257408565831967,
					0.000227525246570192,
					0.426809100739524 ,
					0.685651755353896 ,
					0.177480247177108  ,
					0.0140624536546447 ,
					0.000399662523172413,
					3.53842435837496e-06,
				},
				{
					864.139615974603 ,
					-1107.04797567769 ,
					137.275050612555 ,
					104.753546992766 ,
					7.13892688826536 ,
					0.112263338634835,
					1376.88595178743 ,
					-1776.13593647040 ,
					217.722401334426 ,
					167.369407975834 ,
					11.4154226080259 ,
					0.179609715385726,
					351.985158496693 ,
					-456.670879846764 ,
					55.4637266703440 ,
					42.9096687436616 ,
					2.92973374111487 ,
					0.0461373472193411,
					27.4240022023931 ,
					-35.7733628980449 ,
					4.30853311144907 ,
					3.35409791704109 ,
					0.229313017217905 ,
					0.00361570359281834,
					0.763423292925049 ,
					-1.00111494394894 ,
					0.119629334479884 ,
					0.0937126021760722 ,
					0.00641724028752521,
					0.000101342627733095,
					0.00659948149960499,
					-0.00869943656623366,
					0.00103174926491151,
					0.000813368433858544,
					5.58002443643615e-05 ,
					8.82831342691263e-07,
				},
				{
					2205.71130141894,
					-776.141558841396,
					-1113.28448376848,
					-154.344517152493,
					-6.16604224841230 ,
					-0.0697244317600459,
					-876.497760579035 ,
					307.440602875033 ,
					442.040054083740 ,
					61.3127822286085  ,
					2.45044916002302 ,
					0.0277211363409050,
					-1149.31067038181 ,
					404.058158070479  ,
					580.341481943926  ,
					80.5269955683888 ,
					3.22033280770079 ,
					0.0364575621426572,
					-156.371362914815 ,
					55.0016584266562  ,
					79.0282858822369 ,
					10.9741010056404 ,
					0.439287525442624,
					0.00497880573069235,
					-6.14280520524933 ,
					2.16159686264520 ,
					3.10782106536483  ,
					0.431989704512825 ,
					0.0173140728724888,
					0.000196520520292115 ,
					-0.0682339175385593,
					0.0240220946555860,
					0.0345650504478219 ,
					0.00481028867112797,
					0.000193084897984501,
					2.19537542400338e-06
				}
			};

			for (int i = 0; i < ddR3.GetLength(0); i++)
			{
				for (int j = 0; j < ddR3.GetLength(1); j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedDDR3[i, j], ddR3[i, j], 1e-9));
				}
			}

		}

		[Fact]
		private void TestCollocation2DNaturalSecondDerivatives()
		{
			var nurbs = new NURBS2D(Element.Patch.DegreeKsi, Element.Patch.DegreeHeta,
				Element.Patch.KnotValueVectorKsi, Element.Patch.KnotValueVectorHeta,
				Element.CollocationPoint, Element.ControlPoints);

			var jacobianMatrix = Element.CalculateJacobianMatrix(Element, nurbs);
			var squareDerivatives = Element.CalculateSquareDerivatives(jacobianMatrix);

			var inverseJacobian = jacobianMatrix.Invert();
			var dR = Element.CalculateNaturalDerivatives(nurbs, inverseJacobian);

			var hessianMatrix = Element.CalculateHessian(Element, nurbs, 0);
			var ddR = Element.CalculateNaturalSecondDerivatives(nurbs, hessianMatrix, dR, squareDerivatives);

			var expectedDDR = new double[3, 36]
			{
				{
					88.2335831437404  ,
					-119.505955167848 ,
					18.8302701221316  ,
					12.0693740553013  ,
					0.808164925342631 ,
					0.0126190908235043 ,
					154.724284830929  ,
					-197.204707619928 ,
					22.5376695663091  ,
					18.2677511126776  ,
					1.25175938221591  ,
					0.0197312788429839,
					42.6235614868057  ,
					-51.7028418534275 ,
					4.22675757562360 ,
					4.47257578565579  ,
					0.312813339681897 ,
					0.00497276414998756,
					3.54863847598925 ,
					-4.11226045574429 ,
					0.220907486670591 ,
					0.334576926801204 ,
					0.0238794845927093 ,
					0.000382833440896306, 
					0.105038606513179 ,
					-0.116491200572174 ,
					0.00330592088722368 ,
					0.00895036955848703 ,
					0.000652185207602499 , 
					1.05469277294500e-05 ,
					0.000962133048511820 ,
					-0.00102212901953845 ,
					5.00971639454230e-06 ,
					7.43606226745050e-05  ,
					5.53600016863964e-06 ,
					9.03330287161375e-08,
				},
				{
					1600.90775117971 ,
					2585.93188582269 ,
					671.346629992730  ,
					53.2649005843737  ,
					1.51400418928175 ,
					0.0133939584624100, 
					-2074.81584518481 ,
					-3321.06195621426 ,
					-855.896726173209 ,
					-67.4568119711584  ,
					-1.90539772539023 ,
					-0.0167550678326542,
					257.025103244050  ,
					409.806390644625  ,
					105.298428337080  ,
					8.27762783772177 ,
					0.233270778213823 ,
					0.00204688714618250, 
					192.648911576717  ,
					308.664188519336  ,
					79.6537504329083  ,
					6.28852824369768  ,
					0.177988270646630 ,
					0.00156872535430612,
					12.9388953088588  ,
					20.7501222407440  ,
					5.36100632339471  ,
					0.423850703693734 ,
					0.0120166844216029,
					0.000106108680256955, 
					0.200129129070908 ,
					0.321168095452446 ,
					0.0830637679754102,
					0.00657633399403516,
					0.000186763318487482,
					1.65231149755316e-06,
				},
				{
					-489.494070016594 ,
					194.431782612362 ,
					258.056869285444  ,
					35.5416797178446  ,
					1.41724247413502  ,
					0.0160180548956297,
					185.190946829277  ,
					-94.8844611253714 ,
					-106.542940892258 ,
					-14.3263919541961 ,
					-0.564038770631922 ,
					-0.00631586810305331, 
					263.394903176895 ,
					-89.2000033990269 ,
					-131.167855167361 ,
					-18.2249952952833 ,
					-0.728804373825343 ,
					-0.00824730781331645, 
					36.9595019671664  ,
					-10.3057367660306  ,
					-17.3935757870770  ,
					-2.44802377613172 ,
					-0.0984667985240082,
					-0.00111862268645221, 
					1.49037514272202  ,
					-0.342723546822060 ,
					-0.667995375932686 ,
					-0.0951189302412556 ,
					-0.00384653859600586 ,
					-4.38603822621534e-05 ,
					0.0169599442127733 ,
					-0.00315522958209718, 
					-0.00726085056072592 , 
					-0.00104593923958584 ,
					-4.25259198527444e-05 ,
					-4.86762558078115e-07
				}
			};

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 36; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedDDR[i, j], ddR[i, j], 1e-9));
				}
			}
		}
		
		[Fact]
		private void TestCollocationPointStiffness()
		{
			var stiffnessMatrix = Element.StiffnessMatrix(Element);


			var expectedStiffnessMatrix = new double[2, 72]
			{
				{71269373.1930371,-34963862.1440425,86326396.1395706,13887984.4723116,27890284.6834711,18432633.5203888,3374954.86371781,2538691.40841747,147040.262812224,101231.605295359,1901.86552586240,1144.14677825927,-62797940.7674455,13227924.7735198,-149403999.153288,-6777461.50895510,-30442437.8675070,-7610210.06373269,-587047.590904161,-1023313.71101401,64271.4481680577,-40288.4836165659,1523.84671445658,-451.133435932379,14569488.7496949,18813921.6554925,10080153.2826584,-6371428.81421621,4514418.40589027,-9369132.51195437,809862.146028396,-1301785.37823452,43347.0452809599,-52057.4552732388,625.184027599059,-589.093415236889,7799533.79426816,2639964.42622617,11419802.8050575,-736124.054716472,3087881.33386687,-1242398.27050550,278633.166164329,-174858.841152266,9469.82190319010,-7033.34275171488,102.405199439939,-79.9016204608721,509192.523583931,106455.367337287,785280.393811892,-24480.2533444328,206555.838909382,-47713.9554237633,17285.5072364059,-6794.20930294683,533.848874193793,-274.752756857562,5.24010613399828,-3.13288444729667,7803.00310146480,1211.42458662667,12240.2971855844,-225.373541578370,3195.31082503166,-518.632182908995,261.107419844705,-74.7099456847031,7.79155622409431,-3.03756570376746,0.0734771486659058,-0.0347687541484368},
				{ -34963862.1440425,179317528.052749,13887984.4723116,279571956.210323,18432633.5203888,74498596.1027996,2538691.40841747,6317492.47293727,101231.605295359,197457.353093590,1144.14677825927,1957.21321435567,13227924.7735198,-222050807.197141,-6777461.50895510,-372536659.767168,-7610210.06373269,-93187751.8489012,-1023313.71101401,-6710230.66832101,-40288.4836165659,-161239.773803809,-451.133435932379,-1082.32090523185,18813921.6554925,29883884.5894980,-6371428.81421621,43045098.4610907,-9369132.51195437,11733823.4602800,-1301785.37823452,1081651.57831882,-52057.4552732388,37665.4337475261,-589.093415236889,416.192813041555,2639964.42622617,21306696.1586059,-736124.054716472,33760977.7318489,-1242398.27050550,8761655.83002671,-174858.841152266,703915.403085506,-7033.34275171488,20477.5923356130,-79.9016204608721,187.111764683498,106455.367337287,1425896.57375148,-24480.2533444328,2275752.78247733,-47713.9554237633,589248.724802774,-6794.20930294683,46921.2453889235,-274.752756857562,1345.59881805096,-3.13288444729667,12.0659456002487,1211.42458662667,22029.2171030645,-225.373541578370,35253.8846478690,-518.632182908995,9128.07927210420,-74.7099456847031,725.534089227609,-3.03756570376746,20.7363646754402,-0.0347687541484368,0.185047039297121}
			};

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 72; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedStiffnessMatrix[i, j], stiffnessMatrix[i, j], 1e-9));
				}
			}
		}

		[Fact]
		private void TestCollocationPointCreation()
		{
			var model = new CollocationModel();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\PlateWithHole.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateCollocationModelFromFile();

            //var solverBuilder = new SuiteSparseSolver.Builder();
            //solverBuilder.DofOrderer = new DofOrderer(
            //    new NodeMajorDofOrderingStrategy(), new NullReordering());
            var solverBuilder= new GmresSolver.Builder();
            ISolver solver = new GmresSolver(model,
                new AsymmetricDofOrderer(new RowDofOrderingStrategy()), 
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()));

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model,solver,provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.BuildMatrices();

            var k = solver.LinearSystems[0].Matrix;
            
            Matrix<double> kmatlab=MathNet.Numerics.LinearAlgebra.CreateMatrix.Dense<double>(k.NumRows,k.NumColumns);
            for (int i = 0; i < k.NumRows; i++)
            {
                for (int j = 0; j < k.NumColumns; j++)
                {
                    kmatlab[i, j] = k[i, j];
                }
            }
            MatlabWriter.Write("..\\..\\..\\InputFiles\\KcolMsolve.mat", kmatlab,"Ktotal");
            
        }

        [Fact]
        private void TestCollocationPoint2DBoundaryStiffness()
        {
            var stiffnessMatrix = BoundaryElement.StiffnessMatrix(BoundaryElement);


            var expectedStiffnessMatrix = new double[2, 72]
            {
                {
                    942491.130114883, - 1528622.97410973, - 960229.359363204, 5526.82466713629, 0, 0, 0, 0, 0, 0, 0, 0,
                    1541946.85604645, 604278.375904314, -1534772.04676494, 8833.73948502993, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0,
                    402655.741148221, 795032.315762745, -393396.546787819, 2264.28583707860, 0 ,0, 0, 0, 0, 0, 0, 0,
                    31983.9502496700, 108207.034955708, -30729.3930359958, 176.870208957262, 0 ,0, 0, 0, 0, 0, 0, 0,
                    906.755309707459, 4251.32132908476, -857.631046319342, 4.93629607955652, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0,
                    7.97712650859460, 47.2268812420649, -7.43299716737894, 0.0427823537104383, 0, 0, 0, 0, 0, 0, 0, 0,
                },
                {
                    -1782504.44499018 ,369391.148725693, 5526.82466713629, -336146.346449264, 0, 0, 0, 0, 0, 0, 0, 0,
                    706476.688202230, 523905.386804851, 8833.73948502993, -537275.819700648, 0, 0, 0, 0, 0, 0, 0, 0,
                    927931.867544305, 120303.984740049, 2264.28583707860, -137715.859881830, 0, 0, 0, 0, 0, 0, 0, 0,
                    126273.296906741, 8388.99921530011, 176.870208957262 ,-10757.4019654052, 0, 0, 0, 0, 0, 0, 0, 0,
                    4960.78691125707, 207.172382078212, 4.93629607955652, -300.229877383557, 0, 0, 0, 0, 0, 0, 0, 0,
                    55.1061490101538 ,1.56806701121465, 0.0427823537104383, -2.60206045213942, 0, 0, 0, 0, 0, 0, 0, 0,
                },
            };

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 72; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(-expectedStiffnessMatrix[i, j], stiffnessMatrix[i, j], 1e-9));
                }
            }
        }
    }
}
