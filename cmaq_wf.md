# CMAQ WORKFLOW

1- Herşey WPS ile başlıyor. FNL'i indiriyoruz (ungrib, geogrib, metgrib etc...)
2- run WRF -> `real.exe` && `wrf.exe`
3- WRF çıktıları şu adreste: `mnt/disk2/projects/CityAir/wrf/january`
   Bir anda tüm domainler çalışıyor. Bu yüzden tüm domainler `WRF/month` dizinleri altında.

4- Bundan sonraki adım mcip. Neredeyse hepsi mcip'e bağlı.
   `/mnt/ssd2/APPS/CMAQ/PREP/mcip/scripts`
5- İlk önce `run_mcip_v2.py` kodu ay ay çalıştırılıyor.
   `year, month, day = [2015], [1], list(range(1, 32))`
   `dom_num` çok önemli. değiştirmeyi unutma. Buna göre NCOLS, NROWS da değişecek.
6- mcip çalıştıktan sonra sonuçları `/mnt/disk2/projects/CityAir/mcip` 
   altındaki `domain_size` ile adlandırılmış klasörlere yazıyor. `36km` ve `12km` klasörlerinin altında ay isimleri var. `4km` 'nin altında ise region isimleri var. regionların altında aylar mevcut.

7- `January` olan klasörler, `run_mcip.py` ile oluşturulan günlük `mcip`
   çıktılarıdır. CMAQ'de January kullanılacak. `january_monthly` olanlar ise emisyon hesabında kullanılacak. Dosyalar 2. günden başlıyor ve ayların 1. günleri yok.

8- Ardından orjinal mcip kod (`run_mcip.csh.4km_med`) çalışıyor. Bu monthly çalıştırıyor. Bu scripti de run_mcip.py içerisine gömebiliriz. tek script hem daily hem monthly çalışabilir. Bundaki dosya isimlerinde cityair C ve A küçük. Düzeltilebilir.

9- Sonraki adım `run_icon.4km.csh`. `/mnt/APPS/CMAQ/PREP/icon/scripts`.

10- icon scripti `/mnt/disk2/projects/CityAir/icon` altında ICON* dosyalarını oluşturuyor. En son mediterrenean'lar yapılmış.

11- Next stage is `run_bcon.py`. `/mnt/APPS/CMAQ/PREP/bcon/scripts`.

12- Bu script `/mnt/disk2/projects/CityAir/bcon` altında dosyalar oluşturuyor.
    Dosya isimleri `BCON_v532_CityAir_4km_mediterranean_2015_03*` şeklinde. Fakat, CMAQ bunu 03 değil March istiyor. Bu yüzden aynı klasörde r_renamer.R scripti var. Bu script 03'leri March yapıyor.

13- Sıralama: WRF, mcip, icon, bcon

14- `mnt/ssd2/APPS/programs` 'in altında alper hocanın yazdığı kodlar var.

15- `oceans` kodları. Deniz gridi ise sahil şeridinde dalganın kırıldığı 
    yerleri hesaplıyor. `oceans.py` dosyasını çalıştırıyoruz. Bu surfzone.pickle ve opensea.pickle dosyalarını oluşturuyor. Bu dosyalar var ise ilerliyor. `parameters.yaml` dosyasını kullanıyor.

16- Bu işlemler bitince bir netcdf dosyasına yazdırıyoruz. Örneğin,
    `ocean_4km_cityair_mediterranean.nc` dosyası. 218 ve 219 da pickle dump.
    224. satırdaki kod, aynı klasörde hazır bulunan bizim elle hazırladığımız
    `grid_input_fil_orig` dosyasını `grid_input_fil` dosyasına linkliyor.
    `grid_input_fil_orig` dosyası `parameters.yaml` dosyasının içerisinde.
    Örneğin, `input_grid.txt.4km.mediterranean.cityair` dosyası `input_grid.txt
    dosyasına linkleniyor.

17- `grid_input` dosyaları bir seferlik iş. 4 domain içinde bir seferde
    hazırlanabilir. Fakat mcip'den sonra hazırlanmak zorunda. Bu bilgiler
    mcip'den geliyor. `input_grid.txt` dosya adını kullanmak zorundayız; çünkü fortran bu ismi arıyor. Bu nedenle her seferinde
    `input_grid.txt.4km.mediterrenean.cityair` dosyasını `input_grid.txt`
    olarak linkliyoruz.

18- `grid_input` değerlerini `/mnt/disk2/projects/CityAir/mcip/4km/mediterrenean/january` altındaki `GRIDCRO2D` dosyalarından alıyoruz.

19- `ocean_4km_cityair_mediterranean.nc` dosyasını `/mnt/disk2/projects/CityAir/land` klasörü altına kopyalıyoruz. Ocean kısmı, mcip'e bağlı ama diğerlerinden tamamen bağımsız.

20- Bir sonraki aşama `eproc`. `/mnt/ssd2/APPS/programs/eproc`. Alper hoca hızlı çalışması için `parallel_eproc` yapmış. Her `eprocX` klasörü altında parameters.yaml dosyasını düzenlemek gerekiyor.

21- 3. satırı mediterranean  ve january yapıyoruz. 13. satırda `simulation_day` dediği şey 1. günden başla ve `simulation_ndays` kadar yani 4 gün git demek. `eproc1` klasöründeki scripti ilk defa çalıştırdıktan sonra diğer 6'sına kopyalayıp, `simulation_day` ve `simulation_ndays` parametrelerini ayarlıyoruz el ile. Çünkü ilk `eproc1` `scripts` klasörü altında `Aegean_4_EMEP` adında bir klasör ve bunun altında yaklaşık 8GB'lık pickle dosyaları oluşturuyor. İlk çalıştırmada `eproc1` i çalıştırıp diğerlerine kopyalamak gerekli.

22- Şöyle bir şey yapılabilir. `parallel_eproc` altında her bir domain için bir klasör oluşturulabilir. Örneğin, `parallel_eproc/d04` klasörü oluşturup altında 8 tane eproc klasörü oluşturulabilir. Daha sonra `d05` klasörü oluşturup yine onun altında `eproc` klasörleri oluşturulabilir.

23- `eproc` çıktı dosyaları `/mnt/disk2/projects/CityAiremis/4km/aegean/january` altındaki `CMAQ` klasörü altında oluşuyor. `CMAQ` klasörü altında oluşan dosyaları `january` klasörü altına taşıyıp `CMAQ` klasörünü siliyoruz.

24- `eproc` çıktıları normalde julian date formatında. Alper hoca bunu el ile değiştirip 20150302, 20150303... şeklinde düzeltiyor.

25- 
