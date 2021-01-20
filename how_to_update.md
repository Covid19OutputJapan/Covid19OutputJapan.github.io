## How to update the website

以下、更新日付を2021XXXXとする（ファイル、フォルダ名に使用）

0. `Main_Japan.m`をMATLABで実行して、Figure(.pngファイル)を`/image/2021XXXX`以下に保存
コードやデータは`/_archives/2021XXXX`に保存

1. `/index.md`を`2021XXXX.md`にコピー

- 6行目の
``
permalink: index.html
``
を
``
permalink: 2021XXXX.html
``
に変更

- 10行目の
``
## Updated weekly (Last update on January 20, 2021)
``
を（更新日付に合わせる）
``
## Updated on January 20, 2021
``
に変更

2. `/index.md`を
