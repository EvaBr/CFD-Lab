# Post-processing script generated with paraview
try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

poisuille_1_vtk = LegacyVTKReader( FileNames=['poisuille_1.vtk'] )

RenderView1 = GetRenderView()
CellDatatoPointData1 = CellDatatoPointData()

SetActiveSource(poisuille_1_vtk)
DataRepresentation1 = Show()
DataRepresentation1.Representation = 'Outline'
DataRepresentation1.ScaleFactor = 0.2
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation1.SelectionCellFieldDataArrayName = 'pressure'

RenderView1.CameraFocalPoint = [1.0, 1.0, 1.0]
RenderView1.CameraPosition = [1.0, 1.0, 7.6921304299024635]
RenderView1.CenterOfRotation = [1.0, 1.0, 1.0]

SetActiveSource(CellDatatoPointData1)
DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionPointFieldDataArrayName = 'pressure'
DataRepresentation2.SelectionCellFieldDataArrayName = 'pressure'
DataRepresentation2.Representation = 'Outline'
DataRepresentation2.ScaleFactor = 0.2

DataRepresentation1.Visibility = 0
writer = CreateWriter("poisuille.csv", CellDatatoPointData1)
writer.FieldAssociation = "Points"
writer.UpdatePipeline()
del writer
