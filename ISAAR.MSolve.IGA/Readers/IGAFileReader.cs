using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.IGA.Readers
{
    public class IGAFileReader
    {
	    enum Attributes
	    {
			type, noden, elemn, node, belem, set
	    }

	    enum Types
	    {
			plane, surface
	    }

		public Model Model { get;  }
		public string Filename { get; private set; }

	    public IGAFileReader(Model model, string filename)
	    {
		    Model = model;
		    Filename = filename;
	    }

	    private int numberOfElements;
	    private int controlPointIDcounter=0;
	    private int elementIDCounter = 0;
	    private List<Element> elements = new List<Element>();
	    public void CreateTSplineShellsModelFromFile()
	    {
		    char[] delimeters = { ' ', '=', '\t' };
		    Attributes? name = null;

		    String[] text = System.IO.File.ReadAllLines(Filename);

            Model.PatchesDictionary.Add(0, new Patch { ID = 0 });
		    for (int i = 0; i < text.Length; i++)
		    {
				String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
			    if (line.Length == 0) continue;
			    try
			    {
				    name = (Attributes)Enum.Parse(typeof(Attributes), line[0].ToLower());
			    }
			    catch (Exception exception)
			    {
				    throw new KeyNotFoundException("Variable name " + line[0] + " is not found.");
			    }
				switch (name)
			    {
					case Attributes.type:
						Types type;
						try
						{
							type = (Types)Enum.Parse(typeof(Types), line[1].ToLower());
						}
						catch (Exception exception)
						{
							throw new KeyNotFoundException("Variable name " + line[0] + " is not found.");
						}
						if (type == Types.plane)
							Model.NumberOfDimensions = 2;
						else
							Model.NumberOfDimensions = 3;
						break;
					case Attributes.noden:
						break;
					case Attributes.elemn:
						numberOfElements = Int32.Parse(line[1]);
						break;
					case Attributes.node:
						Model.ControlPointsDictionary.Add(controlPointIDcounter++, new ControlPoint
						{
							X = Double.Parse(line[1], CultureInfo.InvariantCulture),
							Y= Double.Parse(line[2], CultureInfo.InvariantCulture),
							Z = Double.Parse(line[3], CultureInfo.InvariantCulture),
							WeightFactor = Double.Parse(line[4], CultureInfo.InvariantCulture)
						});
						break;
					case Attributes.belem:
						var numberOfElementNodes = Int32.Parse(line[1]);
						var elementDegreeKsi = Int32.Parse(line[2]);
						var elementDegreeHeta = Int32.Parse(line[3]);
						i++;
						line = text[i].Split(delimeters);
						int[] connectivity= new int[numberOfElementNodes];
						for (int j = 0; j < numberOfElementNodes; j++)
							connectivity[j] = Int32.Parse(line[j]);

						Matrix2D extractionOperator = new Matrix2D(numberOfElementNodes,
							(elementDegreeKsi+1)*(elementDegreeHeta+1));
						for (int j = 0; j < numberOfElementNodes; j++)
						{
							line = text[++i].Split(delimeters);
							for (int k = 0; k < (elementDegreeKsi+1)*(elementDegreeHeta+1); k++)
							{
								extractionOperator[j, k] = double.Parse(line[k]);
							}
						}

						if (Model.NumberOfDimensions == 2)
						{
							var element=new Element
							{
								ID = elementIDCounter,
								ElementType = new TSplineElement2D
								{
									DegreeKsi = elementDegreeKsi,
									DegreeHeta = elementDegreeHeta,
									ExtractionOperator = extractionOperator
								}
							};
							foreach (int indexControlPoint in connectivity)
								element.ControlPoints.Add(Model.ControlPointsDictionary[indexControlPoint]);
                            Model.ElementsDictionary.Add(elementIDCounter, element);
                            Model.PatchesDictionary[0].ElementsDictionary.Add(elementIDCounter++, element);
                        }
						else
						{
							var element = new Element
							{
								ID = elementIDCounter,
								ElementType = new TSplineKirchhoffLoveShellElement()
								{
									DegreeKsi = elementDegreeKsi,
									DegreeHeta = elementDegreeHeta,
									ExtractionOperator = extractionOperator
								}
							};
							foreach (int indexControlPoint in connectivity)
								element.ControlPoints.Add(Model.ControlPointsDictionary[indexControlPoint]);
							Model.ElementsDictionary.Add(elementIDCounter,element);
                            Model.PatchesDictionary[0].ElementsDictionary.Add(elementIDCounter++, element);

                        }
						break;
					case Attributes.set:
						break;
			    }
			}

            var a = 0;
        }

    }
}
