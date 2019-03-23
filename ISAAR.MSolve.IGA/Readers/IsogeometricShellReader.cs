using ISAAR.MSolve.IGA.Entities;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.IGA.Readers
{
	public class IsogeometricShellReader
	{
		public Model Model { get; }
		public string Filename { get; private set; }

		enum Attributes
		{
			numberofdimensions,
			numberofpatches, numberofinterfaces,
			patchid, interfaceid,
			degreeksi, degreeheta, degreezeta,
			numberofcpksi, numberofcpheta, numberofcpzeta,
			knotvaluevectorksi, knotvaluevectorheta, knotvaluevectorzeta,
			patchcpid,
			cpcoord, thickness, material, end
		}

		public IsogeometricShellReader(Model model, string filename)
		{
			Model = model;
			Filename = filename;
		}


		private int patchID = -1;
		private int numberOfValues = 0;
		int[] localControlPointIDs;
		int counterElementID = 0;
		private int counterCPID;
		Dictionary<int, int[]> ControlPointIDsDictionary= new Dictionary<int, int[]>();

		public void CreateShellModelFromFile()
		{
			char[] delimeters = { ' ', '=', '\t' };
			IsogeometricShellReader.Attributes? name = null;

			String[] text = System.IO.File.ReadAllLines(Filename);

			for (int i = 0; i < text.Length; i++)
			{
				String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
				if (line.Length == 0)
				{
					continue;
				}
				try
				{
					name = (IsogeometricShellReader.Attributes)Enum.Parse(typeof(IsogeometricShellReader.Attributes), line[0].ToLower());
				}
				catch (Exception exception)
				{
					throw new KeyNotFoundException("Variable name " + line[0] + " is not found.");
				}

				switch (name)
				{
					case IsogeometricShellReader.Attributes.numberofdimensions:
						Model.PatchesDictionary[patchID].NumberOfDimensions = Int32.Parse(line[1]);
						break;
					case IsogeometricShellReader.Attributes.thickness:
						Model.PatchesDictionary[patchID].Thickness = Double.Parse(line[1], CultureInfo.InvariantCulture);
						break;
					case IsogeometricShellReader.Attributes.numberofpatches:
						//Model.NumberOfPatches = Int32.Parse(line[1]);
						break;
					case IsogeometricShellReader.Attributes.material:
						Model.PatchesDictionary[patchID].Material = new ElasticMaterial2D(StressState2D.PlaneStrain) { YoungModulus = Double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture) };
						break;
					case IsogeometricShellReader.Attributes.patchid:
						patchID = Int32.Parse(line[1]);
						Model.PatchesDictionary.Add(patchID,new Patch());
						break;
					case IsogeometricShellReader.Attributes.degreeksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Ksi of a patch must be defined after the patchID");
						Model.PatchesDictionary[patchID].DegreeKsi=Int32.Parse(line[1]);
						break;
					case IsogeometricShellReader.Attributes.degreeheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Heta of a patch must be defined after the patchID");
						Model.PatchesDictionary[patchID].DegreeHeta = Int32.Parse(line[1]);
						break;
					case IsogeometricShellReader.Attributes.numberofcpksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Ksi of a patch must be defined after the patchID");
						Model.PatchesDictionary[patchID].NumberOfControlPointsKsi= Int32.Parse(line[1]);
						break;
					case IsogeometricShellReader.Attributes.numberofcpheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Heta of a patch must be defined after the patchID");
						Model.PatchesDictionary[patchID].NumberOfControlPointsHeta= Int32.Parse(line[1]);
						break;
					case IsogeometricShellReader.Attributes.knotvaluevectorksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Ksi of a patch must be defined after the patchID");
						if (Model.PatchesDictionary[patchID].DegreeKsi == 0 || Model.PatchesDictionary[patchID].NumberOfControlPointsKsi == 0)
							throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
						numberOfValues = Model.PatchesDictionary[patchID].DegreeKsi + Model.PatchesDictionary[patchID].NumberOfControlPointsKsi + 1;
						double[] KnotValueVectorKsi = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
							KnotValueVectorKsi[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						Model.PatchesDictionary[patchID].KnotValueVectorKsi= new Vector(KnotValueVectorKsi);
						break;
					case IsogeometricShellReader.Attributes.knotvaluevectorheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Heta of a patch must be defined after the patchID");
						if (Model.PatchesDictionary[patchID].DegreeHeta == 0 || Model.PatchesDictionary[patchID].NumberOfControlPointsHeta == 0)
							throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
						numberOfValues = Model.PatchesDictionary[patchID].DegreeHeta + Model.PatchesDictionary[patchID].NumberOfControlPointsHeta + 1;
						double[] KnotValueVectorHeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
							KnotValueVectorHeta[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						Model.PatchesDictionary[patchID].KnotValueVectorHeta=new Vector(KnotValueVectorHeta);
						break;
					case IsogeometricShellReader.Attributes.patchcpid:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Control Points ID of a patch must be defined after the patchID");
						int numberOfPatchCP = Model.PatchesDictionary[patchID].NumberOfControlPointsKsi *
						                      Model.PatchesDictionary[patchID].NumberOfControlPointsHeta;
						localControlPointIDs = new int[numberOfPatchCP];
						for (int j = 0; j < numberOfPatchCP; j++)
						{
							localControlPointIDs[j] = Int32.Parse(line[j + 1]);
						}
						ControlPointIDsDictionary.Add(patchID, localControlPointIDs);
						break;
					case IsogeometricShellReader.Attributes.cpcoord:
						var numberOfControlPoints = Int32.Parse(line[1]);
						for (int j = 0; j < numberOfControlPoints; j++)
						{
							i++;
							line = text[i].Split(delimeters);
							int controlPointGlobalID = Int32.Parse(line[0]);
							double x = Double.Parse(line[1], CultureInfo.InvariantCulture);
							double y = Double.Parse(line[2], CultureInfo.InvariantCulture);
							double z = Double.Parse(line[3], CultureInfo.InvariantCulture);
							double w = Double.Parse(line[4], CultureInfo.InvariantCulture);
							ControlPoint controlPoint = new ControlPoint()
							{ ID = controlPointGlobalID, X = x, Y = y, Z = z, WeightFactor = w };
							Model.ControlPointsDictionary.Add(controlPointGlobalID, controlPoint);
						}
						break;
					case IsogeometricShellReader.Attributes.end:
						for (int j = 0; j < ControlPointIDsDictionary[patchID].Length; j++)
							((List<ControlPoint>)Model.PatchesDictionary[patchID].ControlPoints).Add( Model.ControlPointsDictionary[ControlPointIDsDictionary[patchID][j]]);

						Model.PatchesDictionary[patchID].CreateNurbsShell();
						foreach (var element in Model.PatchesDictionary[patchID].Elements)
							Model.ElementsDictionary.Add(counterElementID++, element);
						return;
				}
			}
		}


	}
}
