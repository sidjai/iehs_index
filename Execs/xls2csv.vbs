if WScript.Arguments.Count < 1 Then
    WScript.Echo "Error! Please specify the source path and the destination. Usage: xls2csv SourcePath.xls sheetDestination.csv"
    Wscript.Quit
End If
name = Wscript.Arguments.Item(0)
Dim Xcl
Dim mats

Set Xcl = CreateObject("Excel.Application")
Xcl.DisplayAlerts = False
Set mats = Xcl.Workbooks.Open(name)
mats.Application.DisplayAlerts = False
Dim sh

For sh = 1 To mats.Worksheets.Count
   mats.Worksheets(sh).Activate
   mats.SaveAs Replace(name,".xlsx",sh & ".csv"), 6
Next
mats.Close False
Xcl.Quit