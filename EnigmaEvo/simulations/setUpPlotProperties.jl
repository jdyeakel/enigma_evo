fontSize = 8.97;
plotX_mm = 83#170;
plotSize_mm = (plotX_mm,.75*plotX_mm)
plotSizeInches = plotSize_mm./25.4 #single coloumn size approx 83 mm which is approx 3.2677165354 inches
plotSizePt = 72.27 .* plotSizeInches #latex pt aprox. = 1/72.27 inch
plotSizePt
doublePlotSizePt = plotSizePt.*(1,2)
triplePlotSizePt = plotSizePt.*(1,3)

plotColors = :seaborn_colorblind;
unboundedColor = (:grey,.35);
figurePadding = (0,1,1,0)
legendPatchSize = (20,10)
lineWidth = 1;
plotDataType = "pdf";

CM.set_theme!()
#somehow only worked with absolute paths
CM.update_theme!(fonts = (;regular = "D:\\Dokumente\\Praktikum Merced\\enigma_evo\\EnigmaEvo\\Fonts\\texgyretermes-regular.otf", 
                        italic = "D:\\Dokumente\\Praktikum Merced\\enigma_evo\\EnigmaEvo\\Fonts\\texgyretermes-italic.otf",
                        bold = "D:\\Dokumente\\Praktikum Merced\\enigma_evo\\EnigmaEvo\\Fonts\\texgyretermes-bold.otf",
                        bold_italic = "D:\\Dokumente\\Praktikum Merced\\enigma_evo\\EnigmaEvo\\Fonts\\texgyretermes-bolditalic.otf"))