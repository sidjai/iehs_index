If (WScript.Arguments.Count < 3) Then
    WScript.Echo "Error: Need 3 inputs, output file location, sheet number, sheet"
    Wscript.Quit
End If
xlsName = Wscript.Arguments.Item(0)
sh = Int(Wscript.Arguments.Item(1))
shName = Wscript.Arguments.Item(2)
csvName = Replace(xlsName,".xlsx",".csv")
Dim Xcl
Dim mats
Dim slice

Set Xcl = CreateObject("Excel.Application")
Set fso = CreateObject("Scripting.FileSystemObject")
Set slice = Xcl.Workbooks.Open(csvName)




If (sh = 1) Then
   Set mats = Xcl.Workbooks.Add
   slice.Sheets(1).copy mats.Sheets(1)
   mats.Sheets(2).Delete
   mats.Sheets(2).Delete
   mats.Sheets(2).Delete
Else 
   Set mats = Xcl.Workbooks.Open(xlsName)
   Call mats.Sheets.add(,mats.Sheets(mats.Sheets.Count))
   slice.Sheets(1).copy mats.Sheets(sh)
   mats.Sheets(mats.Sheets.Count).Delete
   
End If
mats.Application.DisplayAlerts = False
'If (mats.Worksheets.Count < sh) Then
'   Set sheet = 
'Else
'   Set sheet = mats.Sheets(sh)
'End If

mats.Sheets(sh).Name = shName

'mats.Sheets(sh).UsedRange.EntireColumn.AutoFit()
mats.SaveAs(xlsName)

mats.Close False
slice.Close False

Xcl.Quit
Set slice = Nothing
Set mats = Nothing
fso.DeleteFile(csvName)