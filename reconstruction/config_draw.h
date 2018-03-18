void ReName(TString &name, TString tag = "", bool justTag = false) {
  name.ReplaceAll("TAG", tag);
  if (justTag) return;
  name.ReplaceAll("*", "x");
  name.ReplaceAll("()", "");
  name.ReplaceAll(":", "_vs_");
  name.ReplaceAll("/", "_div_");
}

TString Shorten(TString branch)
{
  auto array = branch.Tokenize(".");
  auto name = ((TObjString *) array -> At(0)) -> GetString();
  auto count = ((TObjString *) array -> At(1)) -> GetString();

  name.Resize(2);
  return name+count;
}
