# Builds the `cameo_lookup` package data object: one row per CAMEO event
# code, with its Goldstein score, human-readable label, root code, and
# QuadClass/PentaClass recoding. Source values are taken verbatim from the
# project's existing CAMEO/Goldstein reference table (O'Brien 2010 scale;
# quad/penta breakpoints as previously used in this project's GDELT
# pipeline). Re-run this script (and `devtools::document()`) whenever the
# lookup table needs to change; end users should not need to re-run it.

library(dplyr)
library(stringr)
library(tibble)
devtools::load_all(".", quiet = TRUE)

cameo_score <- tibble::tribble(
  ~CAMEOEVENTCODE , ~GoldsteinScore ,
  "01"            ,   0.0           , "010"  ,   0.0 , "011"  ,  -0.1 , "012"  ,  -0.4 , "013"  ,   0.4 ,
  "014"           ,   0.0           , "015"  ,   0.0 , "016"  ,  -2.0 , "017"  ,   0.0 , "018"  ,   3.4 ,
  "019"           ,   3.4           , "02"   ,   3.0 , "020"  ,   3.0 , "021"  ,   3.4 , "0211" ,   3.4 ,
  "0212"          ,   3.4           , "0213" ,   3.4 , "0214" ,   3.4 , "022"  ,   3.2 , "023"  ,   3.4 ,
  "0231"          ,   3.4           , "0232" ,   3.4 , "0233" ,   3.4 , "0234" ,   3.4 , "024"  ,  -0.3 ,
  "0241"          ,  -0.3           , "0242" ,  -0.3 , "0243" ,  -0.3 , "0244" ,  -0.3 , "025"  ,  -0.3 ,
  "0251"          ,  -0.3           , "0252" ,  -0.3 , "0253" ,  -0.3 , "0254" ,  -0.3 , "0255" ,  -0.3 ,
  "0256"          ,  -0.3           , "026"  ,   4.0 , "027"  ,   4.0 , "028"  ,   4.0 , "03"   ,   4.0 ,
  "030"           ,   4.0           , "031"  ,   5.2 , "0311" ,   5.2 , "0312" ,   5.2 , "0313" ,   5.2 ,
  "0314"          ,   5.2           , "032"  ,   4.5 , "033"  ,   5.2 , "0331" ,   5.2 , "0332" ,   5.2 ,
  "0333"          ,   5.2           , "0334" ,   6.0 , "034"  ,   7.0 , "0341" ,   7.0 , "0342" ,   7.0 ,
  "0343"          ,   7.0           , "0344" ,   7.0 , "035"  ,   7.0 , "0351" ,   7.0 , "0352" ,   7.0 ,
  "0353"          ,   7.0           , "0354" ,   7.0 , "0355" ,   7.0 , "0356" ,   7.0 , "036"  ,   4.0 ,
  "037"           ,   5.0           , "038"  ,   7.0 , "039"  ,   5.0 , "04"   ,   1.0 , "040"  ,   1.0 ,
  "041"           ,   1.0           , "042"  ,   1.9 , "043"  ,   2.8 , "044"  ,   2.5 , "045"  ,   5.0 ,
  "046"           ,   7.0           , "05"   ,   3.5 , "050"  ,   3.5 , "051"  ,   3.4 , "052"  ,   3.5 ,
  "053"           ,   3.8           , "054"  ,   6.0 , "055"  ,   7.0 , "056"  ,   7.0 , "057"  ,   8.0 ,
  "06"            ,   6.0           , "060"  ,   6.0 , "061"  ,   6.4 , "062"  ,   7.4 , "063"  ,   7.4 ,
  "064"           ,   7.0           , "07"   ,   7.0 , "070"  ,   7.0 , "071"  ,   7.4 , "072"  ,   8.3 ,
  "073"           ,   7.4           , "074"  ,   8.5 , "075"  ,   7.0 , "08"   ,   5.0 , "080"  ,   5.0 ,
  "081"           ,   5.0           , "0811" ,   5.0 , "0812" ,   5.0 , "0813" ,   5.0 , "0814" ,   5.0 ,
  "082"           ,   5.0           , "083"  ,   5.0 , "0831" ,   5.0 , "0832" ,   5.0 , "0833" ,   5.0 ,
  "0834"          ,   5.0           , "084"  ,   7.0 , "0841" ,   7.0 , "0842" ,   7.0 , "085"  ,   7.0 ,
  "086"           ,   9.0           , "0861" ,   9.0 , "0862" ,   9.0 , "0863" ,   9.0 , "087"  ,   9.0 ,
  "0871"          ,   9.0           , "0872" ,   9.0 , "0873" ,   9.0 , "0874" ,  10.0 , "09"   ,  -2.0 ,
  "090"           ,  -2.0           , "091"  ,  -2.0 , "092"  ,  -2.0 , "093"  ,  -2.0 , "094"  ,  -2.0 ,
  "10"            ,  -5.0           , "100"  ,  -5.0 , "101"  ,  -5.0 , "1011" ,  -5.0 , "1012" ,  -5.0 ,
  "1013"          ,  -5.0           , "1014" ,  -5.0 , "102"  ,  -5.0 , "103"  ,  -5.0 , "1031" ,  -5.0 ,
  "1032"          ,  -5.0           , "1033" ,  -5.0 , "1034" ,  -5.0 , "104"  ,  -5.0 , "1041" ,  -5.0 ,
  "1042"          ,  -5.0           , "1043" ,  -5.0 , "1044" ,  -5.0 , "105"  ,  -5.0 , "1051" ,  -5.0 ,
  "1052"          ,  -5.0           , "1053" ,  -5.0 , "1054" ,  -5.0 , "1055" ,  -5.0 , "1056" ,  -5.0 ,
  "107"           ,  -5.0           , "108"  ,  -5.0 , "11"   ,  -2.0 , "110"  ,  -2.0 , "111"  ,  -2.0 ,
  "112"           ,  -2.0           , "1121" ,  -2.0 , "1122" ,  -2.0 , "1123" ,  -2.0 , "1124" ,  -2.0 ,
  "1125"          ,  -2.0           , "113"  ,  -2.0 , "114"  ,  -2.0 , "115"  ,  -2.0 , "116"  ,  -2.0 ,
  "12"            ,  -4.0           , "120"  ,  -4.0 , "121"  ,  -4.0 , "1211" ,  -4.0 , "1212" ,  -4.0 ,
  "122"           ,  -4.0           , "1221" ,  -4.0 , "1222" ,  -4.0 , "1223" ,  -4.0 , "1224" ,  -4.0 ,
  "123"           ,  -4.0           , "1231" ,  -4.0 , "1232" ,  -4.0 , "1233" ,  -4.0 , "1234" ,  -4.0 ,
  "124"           ,  -4.0           , "1241" ,  -4.0 , "1242" ,  -4.0 , "1243" ,  -4.0 , "1244" ,  -4.0 ,
  "1245"          ,  -4.0           , "1246" ,  -4.0 , "125"  ,  -5.0 , "126"  ,  -5.0 , "127"  ,  -5.0 ,
  "128"           ,  -5.0           , "129"  ,  -5.0 , "13"   ,  -6.0 , "130"  ,  -4.4 , "131"  ,  -5.8 ,
  "1311"          ,  -5.8           , "1312" ,  -5.8 , "1313" ,  -5.8 , "132"  ,  -5.8 , "1321" ,  -5.8 ,
  "1322"          ,  -5.8           , "1323" ,  -5.8 , "1324" ,  -5.8 , "133"  ,  -5.8 , "134"  ,  -5.8 ,
  "135"           ,  -5.8           , "136"  ,  -7.0 , "137"  ,  -7.0 , "138"  ,  -7.0 , "1381" ,  -7.0 ,
  "1382"          ,  -7.0           , "1383" ,  -7.0 , "1384" ,  -7.0 , "1385" ,  -7.0 , "139"  ,  -7.0 ,
  "14"            ,  -6.5           , "140"  ,  -6.5 , "141"  ,  -6.5 , "1411" ,  -6.5 , "1412" ,  -6.5 ,
  "1413"          ,  -6.5           , "1414" ,  -6.5 , "142"  ,  -6.5 , "1421" ,  -6.5 , "1422" ,  -6.5 ,
  "1423"          ,  -6.5           , "1424" ,  -6.5 , "143"  ,  -6.5 , "1431" ,  -6.5 , "1432" ,  -6.5 ,
  "1433"          ,  -6.5           , "1434" ,  -6.5 , "144"  ,  -7.5 , "1441" ,  -7.5 , "1442" ,  -7.5 ,
  "1443"          ,  -7.5           , "1444" ,  -7.5 , "145"  ,  -7.5 , "1451" ,  -7.5 , "1452" ,  -7.5 ,
  "1453"          ,  -7.5           , "1454" ,  -7.5 , "15"   ,  -7.2 , "150"  ,  -7.2 , "151"  ,  -7.2 ,
  "152"           ,  -7.2           , "153"  ,  -7.2 , "154"  ,  -7.2 , "16"   ,  -4.0 , "160"  ,  -4.0 ,
  "161"           ,  -4.0           , "162"  ,  -5.6 , "1621" ,  -5.6 , "1622" ,  -5.6 , "1623" ,  -5.6 ,
  "163"           ,  -8.0           , "164"  ,  -7.0 , "165"  ,  -6.5 , "166"  ,  -7.0 , "1661" ,  -7.0 ,
  "1662"          ,  -7.0           , "1663" ,  -7.0 , "17"   ,  -7.0 , "170"  ,  -7.0 , "171"  ,  -9.2 ,
  "1711"          ,  -9.2           , "1712" ,  -9.2 , "172"  ,  -5.0 , "1721" ,  -5.0 , "1722" ,  -5.0 ,
  "1723"          ,  -5.0           , "1724" ,  -5.0 , "173"  ,  -5.0 , "174"  ,  -5.0 , "175"  ,  -9.0 ,
  "18"            ,  -9.0           , "180"  ,  -9.0 , "181"  ,  -9.0 , "182"  ,  -9.5 , "1821" ,  -9.0 ,
  "1822"          ,  -9.0           , "1823" , -10.0 , "183"  , -10.0 , "1831" , -10.0 , "1832" , -10.0 ,
  "1833"          , -10.0           , "184"  ,  -8.0 , "185"  ,  -8.0 , "186"  , -10.0 , "19"   , -10.0 ,
  "190"           , -10.0           , "191"  ,  -9.5 , "192"  ,  -9.5 , "193"  , -10.0 , "194"  , -10.0 ,
  "195"           , -10.0           , "196"  ,  -9.5 , "20"   , -10.0 , "200"  , -10.0 , "201"  ,  -9.5 ,
  "202"           , -10.0           , "203"  , -10.0 , "204"  , -10.0 , "2041" , -10.0 , "2042" , -10.0
)


cameo_label <- cameo_score %>%
  mutate(
    CAMEOLabel = case_when(
      CAMEOEVENTCODE == "01" ~ "Make Public Statement",
      CAMEOEVENTCODE == "010" ~ "Make statement, not specified below",
      CAMEOEVENTCODE == "011" ~ "Decline comment",
      CAMEOEVENTCODE == "012" ~ "Make pessimistic comment",
      CAMEOEVENTCODE == "013" ~ "Make optimistic comment",
      CAMEOEVENTCODE == "014" ~ "Consider policy option",
      CAMEOEVENTCODE == "015" ~ "Acknowledge or claim responsibility",
      CAMEOEVENTCODE == "016" ~ "Deny responsibility",
      CAMEOEVENTCODE == "017" ~ "Engage in symbolic act",
      CAMEOEVENTCODE == "018" ~ "Make empathetic comment",
      CAMEOEVENTCODE == "019" ~ "Express accord",
      CAMEOEVENTCODE == "02" ~ "Appeal",
      CAMEOEVENTCODE ==
        "020" ~ "Make an appeal or request, not specified below",
      CAMEOEVENTCODE ==
        "021" ~ "Appeal for material cooperation, not specified below",
      CAMEOEVENTCODE == "0211" ~ "Appeal for economic cooperation",
      CAMEOEVENTCODE == "0212" ~ "Appeal for military cooperation",
      CAMEOEVENTCODE == "0213" ~ "Appeal for judicial cooperation",
      CAMEOEVENTCODE == "0214" ~ "Appeal for intelligence",
      CAMEOEVENTCODE == "022" ~ "Appeal for diplomatic cooperation",
      CAMEOEVENTCODE == "023" ~ "Appeal for aid, not specified below",
      CAMEOEVENTCODE == "0231" ~ "Appeal for economic aid",
      CAMEOEVENTCODE == "0232" ~ "Appeal for military aid",
      CAMEOEVENTCODE == "0233" ~ "Appeal for humanitarian aid",
      CAMEOEVENTCODE ==
        "0234" ~ "Appeal for military protection or peacekeeping",
      CAMEOEVENTCODE ==
        "024" ~ "Appeal for political reform, not specified below",
      CAMEOEVENTCODE == "0241" ~ "Appeal for change in leadership",
      CAMEOEVENTCODE == "0242" ~ "Appeal for policy change",
      CAMEOEVENTCODE == "0243" ~ "Appeal for rights",
      CAMEOEVENTCODE == "0244" ~ "Appeal for change in institutions, regime",
      CAMEOEVENTCODE == "025" ~ "Appeal to yield, not specified below",
      CAMEOEVENTCODE ==
        "0251" ~ "Appeal for easing of administrative sanctions",
      CAMEOEVENTCODE == "0252" ~ "Appeal for easing of political dissent",
      CAMEOEVENTCODE == "0253" ~ "Appeal for release of persons or property",
      CAMEOEVENTCODE ==
        "0254" ~ "Appeal for easing of economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "0255" ~ "Appeal for target to allow international involvement (non-mediation)",
      CAMEOEVENTCODE ==
        "0256" ~ "Appeal for de-escalation of military engagement",
      CAMEOEVENTCODE == "026" ~ "Appeal to others to meet or negotiate",
      CAMEOEVENTCODE == "027" ~ "Appeal to others to settle dispute",
      CAMEOEVENTCODE == "028" ~ "Appeal to engage in or accept mediation",
      CAMEOEVENTCODE == "03" ~ "Express Intent to Cooperate",
      CAMEOEVENTCODE ==
        "030" ~ "Express intent to cooperate, not specified below",
      CAMEOEVENTCODE ==
        "031" ~ "Express intent to engage in material cooperation, not specified below",
      CAMEOEVENTCODE == "0311" ~ "Express intent to cooperate economically",
      CAMEOEVENTCODE == "0312" ~ "Express intent to cooperate militarily",
      CAMEOEVENTCODE ==
        "0313" ~ "Express intent to cooperate on judicial matters",
      CAMEOEVENTCODE == "0314" ~ "Express intent to cooperate on intelligence",
      CAMEOEVENTCODE ==
        "032" ~ "Express intent to engage in diplomatic cooperation",
      CAMEOEVENTCODE ==
        "033" ~ "Express intent to provide material aid, not specified below",
      CAMEOEVENTCODE == "0331" ~ "Express intent to provide economic aid",
      CAMEOEVENTCODE == "0332" ~ "Express intent to provide military aid",
      CAMEOEVENTCODE == "0333" ~ "Express intent to provide humanitarian aid",
      CAMEOEVENTCODE ==
        "0334" ~ "Express intent to provide military protection or peacekeeping",
      CAMEOEVENTCODE ==
        "034" ~ "Express intent to institute political reform, not specified below",
      CAMEOEVENTCODE == "0341" ~ "Express intent to change leadership",
      CAMEOEVENTCODE == "0342" ~ "Express intent to change policy",
      CAMEOEVENTCODE == "0343" ~ "Express intent to provide rights",
      CAMEOEVENTCODE ==
        "0344" ~ "Express intent to change institutions, regime",
      CAMEOEVENTCODE == "035" ~ "Express intent to yield, not specified below",
      CAMEOEVENTCODE ==
        "0351" ~ "Express intent to ease administrative sanctions",
      CAMEOEVENTCODE == "0352" ~ "Express intent to ease popular dissent",
      CAMEOEVENTCODE ==
        "0353" ~ "Express intent to release persons or property",
      CAMEOEVENTCODE ==
        "0354" ~ "Express intent to ease economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "0355" ~ "Express intent to allow international involvement (non-mediation)",
      CAMEOEVENTCODE ==
        "0356" ~ "Express intent to de-escalate military engagement",
      CAMEOEVENTCODE == "036" ~ "Express intent to meet or negotiate",
      CAMEOEVENTCODE == "037" ~ "Express intent to settle dispute",
      CAMEOEVENTCODE == "038" ~ "Express intent to accept mediation",
      CAMEOEVENTCODE == "039" ~ "Express intent to mediate",
      CAMEOEVENTCODE == "04" ~ "Consult",
      CAMEOEVENTCODE == "040" ~ "Consult, not specified below",
      CAMEOEVENTCODE == "041" ~ "Discuss by telephone",
      CAMEOEVENTCODE == "042" ~ "Make a visit",
      CAMEOEVENTCODE == "043" ~ "Host a visit",
      CAMEOEVENTCODE == "044" ~ "Meet at a third location",
      CAMEOEVENTCODE == "045" ~ "Mediate",
      CAMEOEVENTCODE == "046" ~ "Engage in negotiation",
      CAMEOEVENTCODE == "05" ~ "Engage in Diplomatic Cooperation",
      CAMEOEVENTCODE ==
        "050" ~ "Engage in diplomatic cooperation, not specified below",
      CAMEOEVENTCODE == "051" ~ "Praise or endorse",
      CAMEOEVENTCODE == "052" ~ "Defend verbally",
      CAMEOEVENTCODE == "053" ~ "Rally support on behalf of",
      CAMEOEVENTCODE == "054" ~ "Grant diplomatic recognition",
      CAMEOEVENTCODE == "055" ~ "Apologize",
      CAMEOEVENTCODE == "056" ~ "Forgive",
      CAMEOEVENTCODE == "057" ~ "Sign formal agreement",
      CAMEOEVENTCODE == "06" ~ "Engage in Material Cooperation",
      CAMEOEVENTCODE ==
        "060" ~ "Engage in material cooperation, not specified below",
      CAMEOEVENTCODE == "061" ~ "Cooperate economically",
      CAMEOEVENTCODE == "062" ~ "Cooperate militarily",
      CAMEOEVENTCODE == "063" ~ "Engage in judicial cooperation",
      CAMEOEVENTCODE == "064" ~ "Share intelligence or information",
      CAMEOEVENTCODE == "07" ~ "Provide Aid",
      CAMEOEVENTCODE == "070" ~ "Provide aid, not specified below",
      CAMEOEVENTCODE == "071" ~ "Provide economic aid",
      CAMEOEVENTCODE == "072" ~ "Provide military aid",
      CAMEOEVENTCODE == "073" ~ "Provide humanitarian aid",
      CAMEOEVENTCODE == "074" ~ "Provide military protection or peacekeeping",
      CAMEOEVENTCODE == "075" ~ "Grant asylum",
      CAMEOEVENTCODE == "08" ~ "Yield",
      CAMEOEVENTCODE == "080" ~ "Yield, not specified below",
      CAMEOEVENTCODE ==
        "081" ~ "Ease administrative sanctions, not specified below",
      CAMEOEVENTCODE == "0811" ~ "Ease restrictions on political freedoms",
      CAMEOEVENTCODE == "0812" ~ "Ease ban on political parties or politicians",
      CAMEOEVENTCODE == "0813" ~ "Ease curfew",
      CAMEOEVENTCODE == "0814" ~ "Ease state of emergency or martial law",
      CAMEOEVENTCODE == "082" ~ "Ease political dissent",
      CAMEOEVENTCODE ==
        "083" ~ "Accede to requests or demands for political reform, not specified below",
      CAMEOEVENTCODE == "0831" ~ "Accede to demands for change in leadership",
      CAMEOEVENTCODE == "0832" ~ "Accede to demands for change in policy",
      CAMEOEVENTCODE == "0833" ~ "Accede to demands for rights",
      CAMEOEVENTCODE ==
        "0834" ~ "Accede to demands for change in institutions, regime",
      CAMEOEVENTCODE == "084" ~ "Return, release, not specified below",
      CAMEOEVENTCODE == "0841" ~ "Return, release person(s)",
      CAMEOEVENTCODE == "0842" ~ "Return, release property",
      CAMEOEVENTCODE == "085" ~ "Ease economic sanctions, boycott, embargo",
      CAMEOEVENTCODE ==
        "086" ~ "Allow international involvement, not specified below",
      CAMEOEVENTCODE == "0861" ~ "Receive deployment of peacekeepers",
      CAMEOEVENTCODE == "0862" ~ "Receive inspectors",
      CAMEOEVENTCODE == "0863" ~ "Allow humanitarian access",
      CAMEOEVENTCODE == "087" ~ "De-escalate military engagement",
      CAMEOEVENTCODE == "0871" ~ "Declare truce, ceasefire",
      CAMEOEVENTCODE == "0872" ~ "Ease military blockade",
      CAMEOEVENTCODE == "0873" ~ "Demobilize armed forces",
      CAMEOEVENTCODE == "0874" ~ "Retreat or surrender militarily",
      CAMEOEVENTCODE == "09" ~ "Investigate",
      CAMEOEVENTCODE == "090" ~ "Investigate, not specified below",
      CAMEOEVENTCODE == "091" ~ "Investigate crime, corruption",
      CAMEOEVENTCODE == "092" ~ "Investigate human rights abuses",
      CAMEOEVENTCODE == "093" ~ "Investigate military action",
      CAMEOEVENTCODE == "094" ~ "Investigate war crimes",
      CAMEOEVENTCODE == "10" ~ "Demand",
      CAMEOEVENTCODE == "100" ~ "Demand, not specified below",
      CAMEOEVENTCODE ==
        "101" ~ "Demand material cooperation, not specified below",
      CAMEOEVENTCODE == "1011" ~ "Demand economic cooperation",
      CAMEOEVENTCODE == "1012" ~ "Demand military cooperation",
      CAMEOEVENTCODE == "1013" ~ "Demand judicial cooperation",
      CAMEOEVENTCODE == "1014" ~ "Demand intelligence cooperation",
      CAMEOEVENTCODE == "102" ~ "Demand diplomatic cooperation",
      CAMEOEVENTCODE == "103" ~ "Demand material aid, not specified below",
      CAMEOEVENTCODE == "1031" ~ "Demand economic aid",
      CAMEOEVENTCODE == "1032" ~ "Demand military aid",
      CAMEOEVENTCODE == "1033" ~ "Demand humanitarian aid",
      CAMEOEVENTCODE == "1034" ~ "Demand military protection or peacekeeping",
      CAMEOEVENTCODE == "104" ~ "Demand political reform, not specified below",
      CAMEOEVENTCODE == "1041" ~ "Demand change in leadership",
      CAMEOEVENTCODE == "1042" ~ "Demand policy change",
      CAMEOEVENTCODE == "1043" ~ "Demand rights",
      CAMEOEVENTCODE == "1044" ~ "Demand change in institutions, regime",
      CAMEOEVENTCODE ==
        "105" ~ "Demand that target yields, not specified below",
      CAMEOEVENTCODE == "1051" ~ "Demand easing of administrative sanctions",
      CAMEOEVENTCODE == "1052" ~ "Demand easing of political dissent",
      CAMEOEVENTCODE == "1053" ~ "Demand release of persons or property",
      CAMEOEVENTCODE ==
        "1054" ~ "Demand easing of economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "1055" ~ "Demand that target allows international involvement (non-mediation)",
      CAMEOEVENTCODE == "1056" ~ "Demand de-escalation of military engagement",
      CAMEOEVENTCODE == "107" ~ "Demand settling of dispute",
      CAMEOEVENTCODE == "108" ~ "Demand mediation",
      CAMEOEVENTCODE == "11" ~ "Disapprove",
      CAMEOEVENTCODE == "110" ~ "Disapprove, not specified below",
      CAMEOEVENTCODE == "111" ~ "Criticize or denounce",
      CAMEOEVENTCODE == "112" ~ "Accuse, not specified below",
      CAMEOEVENTCODE == "1121" ~ "Accuse of crime, corruption",
      CAMEOEVENTCODE == "1122" ~ "Accuse of human rights abuses",
      CAMEOEVENTCODE == "1123" ~ "Accuse of aggression",
      CAMEOEVENTCODE == "1124" ~ "Accuse of war crimes",
      CAMEOEVENTCODE == "1125" ~ "Accuse of espionage, treason",
      CAMEOEVENTCODE == "113" ~ "Rally opposition against",
      CAMEOEVENTCODE == "114" ~ "Complain officially",
      CAMEOEVENTCODE == "115" ~ "Bring lawsuit against",
      CAMEOEVENTCODE == "116" ~ "Find guilty or liable (legally)",
      CAMEOEVENTCODE == "12" ~ "Reject",
      CAMEOEVENTCODE == "120" ~ "Reject, not specified below",
      CAMEOEVENTCODE == "121" ~ "Reject material cooperation",
      CAMEOEVENTCODE == "1211" ~ "Reject economic cooperation",
      CAMEOEVENTCODE == "1212" ~ "Reject military cooperation",
      CAMEOEVENTCODE ==
        "122" ~ "Reject request or demand for material aid, not specified below",
      CAMEOEVENTCODE == "1221" ~ "Reject request for economic aid",
      CAMEOEVENTCODE == "1222" ~ "Reject request for military aid",
      CAMEOEVENTCODE == "1223" ~ "Reject request for humanitarian aid",
      CAMEOEVENTCODE ==
        "1224" ~ "Reject request for military protection or peacekeeping",
      CAMEOEVENTCODE ==
        "123" ~ "Reject request or demand for political reform, not specified below",
      CAMEOEVENTCODE == "1231" ~ "Reject request for change in leadership",
      CAMEOEVENTCODE == "1232" ~ "Reject request for policy change",
      CAMEOEVENTCODE == "1233" ~ "Reject request for rights",
      CAMEOEVENTCODE ==
        "1234" ~ "Reject request for change in institutions, regime",
      CAMEOEVENTCODE == "124" ~ "Refuse to yield, not specified below",
      CAMEOEVENTCODE == "1241" ~ "Refuse to ease administrative sanctions",
      CAMEOEVENTCODE == "1242" ~ "Refuse to ease popular dissent",
      CAMEOEVENTCODE == "1243" ~ "Refuse to release persons or property",
      CAMEOEVENTCODE ==
        "1244" ~ "Refuse to ease economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "1245" ~ "Refuse to allow international involvement (non-mediation)",
      CAMEOEVENTCODE == "1246" ~ "Refuse to de-escalate military engagement",
      CAMEOEVENTCODE ==
        "125" ~ "Reject proposal to meet, discuss, or negotiate",
      CAMEOEVENTCODE == "126" ~ "Reject mediation",
      CAMEOEVENTCODE == "127" ~ "Reject plan, agreement to settle dispute",
      CAMEOEVENTCODE == "128" ~ "Defy norms, law",
      CAMEOEVENTCODE == "129" ~ "Veto",
      CAMEOEVENTCODE == "13" ~ "Threaten",
      CAMEOEVENTCODE == "130" ~ "Threaten, not specified below",
      CAMEOEVENTCODE == "131" ~ "Threaten non-force, not specified below",
      CAMEOEVENTCODE == "1311" ~ "Threaten to reduce or stop aid",
      CAMEOEVENTCODE == "1312" ~ "Threaten with sanctions, boycott, embargo",
      CAMEOEVENTCODE == "1313" ~ "Threaten to reduce or break relations",
      CAMEOEVENTCODE ==
        "132" ~ "Threaten with administrative sanctions, not specified below",
      CAMEOEVENTCODE ==
        "1321" ~ "Threaten with restrictions on political freedoms",
      CAMEOEVENTCODE ==
        "1322" ~ "Threaten to ban political parties or politicians",
      CAMEOEVENTCODE == "1323" ~ "Threaten to impose curfew",
      CAMEOEVENTCODE ==
        "1324" ~ "Threaten to impose state of emergency or martial law",
      CAMEOEVENTCODE == "133" ~ "Threaten with political dissent, protest",
      CAMEOEVENTCODE == "134" ~ "Threaten to halt negotiations",
      CAMEOEVENTCODE == "135" ~ "Threaten to halt mediation",
      CAMEOEVENTCODE ==
        "136" ~ "Threaten to halt international involvement (non-mediation)",
      CAMEOEVENTCODE == "137" ~ "Threaten with repression",
      CAMEOEVENTCODE ==
        "138" ~ "Threaten with military force, not specified below",
      CAMEOEVENTCODE == "1381" ~ "Threaten blockade",
      CAMEOEVENTCODE == "1382" ~ "Threaten occupation",
      CAMEOEVENTCODE == "1383" ~ "Threaten unconventional violence",
      CAMEOEVENTCODE == "1384" ~ "Threaten conventional attack",
      CAMEOEVENTCODE == "1385" ~ "Threaten attack with WMD",
      CAMEOEVENTCODE == "139" ~ "Give ultimatum",
      CAMEOEVENTCODE == "14" ~ "Protest",
      CAMEOEVENTCODE ==
        "140" ~ "Engage in political dissent, not specified below",
      CAMEOEVENTCODE == "141" ~ "Demonstrate or rally, not specified below",
      CAMEOEVENTCODE == "1411" ~ "Demonstrate for leadership change",
      CAMEOEVENTCODE == "1412" ~ "Demonstrate for policy change",
      CAMEOEVENTCODE == "1413" ~ "Demonstrate for rights",
      CAMEOEVENTCODE ==
        "1414" ~ "Demonstrate for change in institutions, regime",
      CAMEOEVENTCODE == "142" ~ "Conduct hunger strike, not specified below",
      CAMEOEVENTCODE == "1421" ~ "Conduct hunger strike for leadership change",
      CAMEOEVENTCODE == "1422" ~ "Conduct hunger strike for policy change",
      CAMEOEVENTCODE == "1423" ~ "Conduct hunger strike for rights",
      CAMEOEVENTCODE ==
        "1424" ~ "Conduct hunger strike for change in institutions, regime",
      CAMEOEVENTCODE ==
        "143" ~ "Conduct strike or boycott, not specified below",
      CAMEOEVENTCODE ==
        "1431" ~ "Conduct strike or boycott for leadership change",
      CAMEOEVENTCODE == "1432" ~ "Conduct strike or boycott for policy change",
      CAMEOEVENTCODE == "1433" ~ "Conduct strike or boycott for rights",
      CAMEOEVENTCODE ==
        "1434" ~ "Conduct strike or boycott for change in institutions, regime",
      CAMEOEVENTCODE == "144" ~ "Obstruct passage, block, not specified below",
      CAMEOEVENTCODE == "1441" ~ "Obstruct passage to demand leadership change",
      CAMEOEVENTCODE == "1442" ~ "Obstruct passage to demand policy change",
      CAMEOEVENTCODE == "1443" ~ "Obstruct passage to demand rights",
      CAMEOEVENTCODE ==
        "1444" ~ "Obstruct passage to demand change in institutions, regime",
      CAMEOEVENTCODE == "145" ~ "Protest violently, riot, not specified below",
      CAMEOEVENTCODE ==
        "1451" ~ "Engage in violent protest for leadership change",
      CAMEOEVENTCODE == "1452" ~ "Engage in violent protest for policy change",
      CAMEOEVENTCODE == "1453" ~ "Engage in violent protest for rights",
      CAMEOEVENTCODE ==
        "1454" ~ "Engage in violent protest for change in institutions, regime",
      CAMEOEVENTCODE == "15" ~ "Exhibit Force Posture",
      CAMEOEVENTCODE ==
        "150" ~ "Demonstrate military or police power, not specified below",
      CAMEOEVENTCODE == "151" ~ "Increase police alert status",
      CAMEOEVENTCODE == "152" ~ "Increase military alert status",
      CAMEOEVENTCODE == "153" ~ "Mobilize or increase police power",
      CAMEOEVENTCODE == "154" ~ "Mobilize or increase armed forces",
      CAMEOEVENTCODE == "155" ~ "Mobilize or increase cyber-forces",
      CAMEOEVENTCODE == "16" ~ "Reduce Relations",
      CAMEOEVENTCODE == "160" ~ "Reduce relations, not specified below",
      CAMEOEVENTCODE == "161" ~ "Reduce or break diplomatic relations",
      CAMEOEVENTCODE ==
        "162" ~ "Reduce or stop material aid, not specified below",
      CAMEOEVENTCODE == "1621" ~ "Reduce or stop economic assistance",
      CAMEOEVENTCODE == "1622" ~ "Reduce or stop military assistance",
      CAMEOEVENTCODE == "1623" ~ "Reduce or stop humanitarian assistance",
      CAMEOEVENTCODE == "163" ~ "Impose embargo, boycott, or sanctions",
      CAMEOEVENTCODE == "164" ~ "Halt negotiations",
      CAMEOEVENTCODE == "165" ~ "Halt mediation",
      CAMEOEVENTCODE == "166" ~ "Expel or withdraw, not specified below",
      CAMEOEVENTCODE == "1661" ~ "Expel or withdraw peacekeepers",
      CAMEOEVENTCODE == "1662" ~ "Expel or withdraw inspectors, observers",
      CAMEOEVENTCODE == "1663" ~ "Expel or withdraw aid agencies",
      CAMEOEVENTCODE == "17" ~ "Coerce",
      CAMEOEVENTCODE == "170" ~ "Coerce, not specified below",
      CAMEOEVENTCODE == "171" ~ "Seize or damage property, not specified below",
      CAMEOEVENTCODE == "1711" ~ "Confiscate property",
      CAMEOEVENTCODE == "1712" ~ "Destroy property",
      CAMEOEVENTCODE ==
        "172" ~ "Impose administrative sanctions, not specified below",
      CAMEOEVENTCODE == "1721" ~ "Impose restrictions on political freedoms",
      CAMEOEVENTCODE == "1722" ~ "Ban political parties or politicians",
      CAMEOEVENTCODE == "1723" ~ "Impose curfew",
      CAMEOEVENTCODE == "1724" ~ "Impose state of emergency or martial law",
      CAMEOEVENTCODE == "173" ~ "Arrest, detain, or charge with legal action",
      CAMEOEVENTCODE == "174" ~ "Expel or deport individuals",
      CAMEOEVENTCODE == "175" ~ "Use tactics of violent repression",
      CAMEOEVENTCODE == "176" ~ "Attack cybernetically",
      CAMEOEVENTCODE == "18" ~ "Assault",
      CAMEOEVENTCODE ==
        "180" ~ "Use unconventional violence, not specified below",
      CAMEOEVENTCODE == "181" ~ "Abduct, hijack, or take hostage",
      CAMEOEVENTCODE == "182" ~ "Physically assault, not specified below",
      CAMEOEVENTCODE == "1821" ~ "Sexually assault",
      CAMEOEVENTCODE == "1822" ~ "Torture",
      CAMEOEVENTCODE == "1823" ~ "Kill by physical assault",
      CAMEOEVENTCODE ==
        "183" ~ "Conduct suicide, car, or other non-military bombing, not specified below",
      CAMEOEVENTCODE == "1831" ~ "Carry out suicide bombing",
      CAMEOEVENTCODE == "1832" ~ "Carry out vehicular bombing",
      CAMEOEVENTCODE == "1833" ~ "Carry out roadside bombing",
      CAMEOEVENTCODE == "1834" ~ "Carry out location bombing",
      CAMEOEVENTCODE == "184" ~ "Use as human shield",
      CAMEOEVENTCODE == "185" ~ "Attempt to assassinate",
      CAMEOEVENTCODE == "186" ~ "Assassinate",
      CAMEOEVENTCODE == "19" ~ "Fight",
      CAMEOEVENTCODE ==
        "190" ~ "Use conventional military force, not specified below",
      CAMEOEVENTCODE == "191" ~ "Impose blockade, restrict movement",
      CAMEOEVENTCODE == "192" ~ "Occupy territory",
      CAMEOEVENTCODE == "193" ~ "Fight with small arms and light weapons",
      CAMEOEVENTCODE == "194" ~ "Fight with artillery and tanks",
      CAMEOEVENTCODE == "195" ~ "Employ aerial weapons, not specified below",
      CAMEOEVENTCODE == "1951" ~ "Employ precision-guided aerial munitions",
      CAMEOEVENTCODE == "1952" ~ "Employ remotely piloted aerial munitions",
      CAMEOEVENTCODE == "196" ~ "Violate ceasefire",
      CAMEOEVENTCODE == "20" ~ "Use Unconventional Mass Violence",
      CAMEOEVENTCODE ==
        "200" ~ "Use unconventional mass violence, not specified below",
      CAMEOEVENTCODE == "201" ~ "Engage in mass expulsion",
      CAMEOEVENTCODE == "202" ~ "Engage in mass killings",
      CAMEOEVENTCODE == "203" ~ "Engage in ethnic cleansing",
      CAMEOEVENTCODE ==
        "204" ~ "Use weapons of mass destruction, not specified below",
      CAMEOEVENTCODE ==
        "2041" ~ "Use chemical, biological, or radiological weapons",
      CAMEOEVENTCODE == "2042" ~ "Detonate nuclear weapons",
      TRUE ~ NA_character_
    )
  )

cameo_lookup <- cameo_label %>%
  mutate(
    root = get_root(CAMEOEVENTCODE),
    QuadClass = assign_quad(root),
    PentaClass = assign_penta(root, QuadClass),
    PentaClass_modified = if_else(
      GoldsteinScore <= 1 & QuadClass == 1,
      0L,
      PentaClass
    )
  ) %>%
  select(-root) %>%
  select(CAMEOEVENTCODE, CAMEOLabel, everything()) %>%
  distinct(CAMEOEVENTCODE, .keep_all = TRUE)

usethis::use_data(cameo_lookup, overwrite = TRUE)
