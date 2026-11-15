# راهنمای انتشار PlasmaROM-H2 در GitHub

پروژه در ZIP آماده انتشار است، اما **تاکنون به حساب GitHub شما Push نشده است**.

1. در GitHub یک مخزن جدید با نام `PlasmaROM-H2` بسازید؛ بهتر است آن را Public بگذارید. هنگام ساخت، گزینه‌های خودکار README، License و .gitignore را انتخاب نکنید، چون داخل پروژه وجود دارند.
2. ZIP را استخراج کنید و در ترمینال وارد پوشه `PlasmaROM-H2` شوید.
3. دستورات زیر را اجرا کنید و `YOUR_USERNAME` را با نام کاربری واقعی خود جایگزین کنید:

```bash
git init
git branch -M main
git add .
git commit -m "Initial release: PlasmaROM-H2 v0.12.0"
git remote add origin https://github.com/YOUR_USERNAME/PlasmaROM-H2.git
git push -u origin main
```

اگر Git نام و ایمیل نداشت، قبل از `git commit` از `git config --global user.name` و `git config --global user.email` برای تنظیم هویت خود استفاده کنید. هنگام Push، از راه‌های احراز هویت GitHub مانند GitHub CLI یا Credential Manager استفاده کنید؛ رمز عبور GitHub را داخل فایل‌ها ذخیره نکنید.

4. پس از Push در برگه **Actions** ابتدا نتیجه CI عددی را ببینید؛ سپس Workflow شیمی واقعی را جدا بررسی کنید. اگر Cantera یا مکانیزم‌ها نصب/یافت نشوند، Workflow باید شکست را ثبت کند و نباید ادعای موفقیت فیزیکی کرد.
5. اگر می‌خواهید انتشار نسخه‌ای بسازید، در Releases یک تگ `v0.12.0` و عنوان `PlasmaROM-H2 v0.12.0 — Numerical prototype` قرار دهید. در توضیح انتشار اشاره کنید «real PAC and experimental validation not completed».

**توجه:** عبور از تست‌های عددی یا یکپارچه Cantera، به‌تنهایی اعتبارسنجی تجربی Plasma-Assisted Combustion نیست. محدودیت‌ها در `docs/SCIENTIFIC_SCOPE.md` آمده‌اند.
