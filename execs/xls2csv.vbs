if WScript.Arguments.Count < 2 Then
    WScript.Echo "Error! Please specify the source path and the allSheetsFlag"
    Wscript.Quit
End If
name = Wscript.Arguments.Item(0)
inMax = Int(Wscript.Arguments.Item(1))

Dim Xcl
Dim mats

Set Xcl = CreateObject("Excel.Application")
Xcl.DisplayAlerts = False
Set mats = Xcl.Workbooks.Open(name)
mats.Application.DisplayAlerts = False
Dim sh

realMax = mats.Worksheets.Count
if(inMax > realMax Or inMax = 9999) Then
   parseMax = realMax
Else 
   parseMax = inMax
End If

For sh = 1 To parseMax 
   mats.Worksheets(sh).Activate
   mats.SaveAs Replace(name,".xlsx",sh & ".csv"), 6
Next
mats.Close False
Xcl.Quit