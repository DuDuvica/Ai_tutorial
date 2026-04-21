void collinsSoperFrame() {
  TCanvas *c = new TCanvas("c","CS frame",800,800);
  c->Range(-1.8,-1.4,1.9,2.0);

  TArrow *a1 = new TArrow(0,0,-0.9,1.2,0.03,"|>");
  TArrow *a2 = new TArrow(0,0, 0.9,1.2,0.03,"|>");
  TArrow *az = new TArrow(0,0, 0.0,1.5,0.03,"|>");
  TArrow *ax = new TArrow(0,0, 1.3,0.0,0.03,"|>");
  TArrow *al = new TArrow(0,0, 0.9,1.0,0.03,"|>");
  TArrow *alp = new TArrow(0,0,-0.9,-1.0,0.03,"|>");

  a1->SetLineColor(kBlue+1); a2->SetLineColor(kBlue+1);
  al->SetLineColor(kRed+1);
  alp->SetLineColor(kOrange+7);

  a1->Draw(); a2->Draw(); az->Draw(); ax->Draw(); al->Draw(); alp->Draw();

  TLatex t; t.SetTextSize(0.03);
  t.DrawLatex(-1.55,1.8,"Collins-Soper frame");
  t.DrawLatex(-1.0,1.25,"p_{1}");
  t.DrawLatex(0.95,1.25,"p_{2}");
  t.DrawLatex(0.05,1.52,"z_{CS}");
  t.DrawLatex(1.32,-0.02,"x_{CS}");
  t.DrawLatex(0.92,1.03,"l^{-}");
  t.DrawLatex(-1.03,-1.05,"l^{+}");

  TEllipse *theta = new TEllipse(0,0,0.35,0.35,48,90);
  theta->SetFillStyle(0); theta->SetLineColor(kGreen+2); theta->Draw();
  t.SetTextColor(kGreen+2); t.DrawLatex(0.12,0.42,"#theta_{CS}");

  TEllipse *phi = new TEllipse(0,0,0.35,0.35,0,48.2);
  phi->SetFillStyle(0); phi->SetLineStyle(2); phi->SetLineColor(kMagenta+1); phi->Draw();
  t.SetTextColor(kMagenta+1); t.DrawLatex(0.45,0.08,"#phi_{CS}");
}